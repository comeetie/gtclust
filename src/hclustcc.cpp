#include <Rcpp.h>
#include "GTMethod.h"
#include "node.h"
#include "LaplaceDirichlet.h"
#include "cholmod_tools.h"
using namespace Rcpp;




double cut_cost(cholmod_factor * Lintra,bayesian_node * node_g,bayesian_node * node_h,int * permutation,cholmod_common * com)
{
  
  std::vector<std::pair<int, int>> cutset = node_g->neibs.at(node_h->id).edges;
  double logdetdiff = 0;
  if(cutset.size()>1){
    int pivot = node_g->intra_pivot; 
    int old_pivot=node_h->intra_pivot;
    std::set<int> old_pivot_edges = node_h->intra_pivot_edges;
    double logdet = cholmod_tools_cutcost(Lintra,pivot,old_pivot,old_pivot_edges,&cutset,permutation,&node_g->intra_nodes,&node_h->intra_nodes,com);
    logdetdiff = logdet-node_g->lognbtree-node_h->lognbtree-log(cutset.size());
    //Rcout << "logdet :: " << logdet << std::endl;
    //Rcout << "logdetdiff :: " << logdetdiff << std::endl;
  }
  
  return logdetdiff;
}



std::map<int, int> build_inter_links(bayesian_node * node,const std::vector<bayesian_node> * graph){
  
  //Rcout << "Node : " << node->cl << std::endl;
  int diag_val=0;
  std::map<int, int> inter_links;
  for (auto itr = node->neibs.begin(); itr!=node->neibs.end();++itr){
    int clto = graph->at(itr->first).i_inter;
    multiedge e = itr->second;
    inter_links.insert(std::make_pair(clto,-e.size));
    diag_val+=e.size;
  }
  inter_links.insert(std::make_pair(node->i_inter,diag_val));
  return(inter_links);
}

// void print_inter(std::map<int, int> inter){
//   
//   for (auto itr = inter.begin(); itr!=inter.end();++itr){
//     Rcout << "(" << itr->first << " , " << itr->second << "), ";
//   }
//   Rcout << std::endl;
// }
// 
// 
// 
// void print_pq(std::multimap<double,std::pair<int, int>,std::less<double>> priority_queue){
//   for(auto it=priority_queue.begin();it!=priority_queue.end();it++){
//     std::pair<int,int> edge = it->second;
//     Rcout << edge.first << " - " << edge.second << ":" << it->first << std::endl;
//   }
// }

// void print_pqhead(std::multimap<double,std::pair<int, int>,std::less<double>> priority_queue){
//   int i =0;
//   for(auto it=priority_queue.begin();it!=priority_queue.end();it++){
//     std::pair<int,int> edge = it->second;
//     Rcout << edge.first << " - " << edge.second << ":" << it->first << std::endl;
//     if(i>10){
//       break;
//     }
//     i++;
//   }
// }


GTMethod::GTMethod * init_method(List method_obj){
  GTMethod::GTMethod * method;
  if(!method_obj.inherits("gtmethod")){
    stop("Method should be a gtmethod object.");
  }
  std::string method_name = method_obj["method"];
  if(method_name=="ward"){
    method = new GTMethod::ward();
  }else if(method_name=="centroid"){
    method = new GTMethod::centroid();
  }else if(method_name=="median"){
    method = new GTMethod::median();
  }else if(method_name=="chisq"){
    method = new GTMethod::chisq(); 
  }else if(method_name=="bayes_mom"){
    double beta = method_obj["beta"];
    method = new GTMethod::bayes_mom(beta); 
  }else if(method_name=="bayes_dgmm"){
    double kappa = method_obj["kappa"];
    double tau = method_obj["tau"];
    double beta = method_obj["beta"];
    NumericVector mu = method_obj["mu"];
    method = new GTMethod::bayes_dgmm(kappa,tau,beta,mu); 
  }else if(method_name=="bayes_dirichlet"){
    NumericVector lambda = method_obj["lambda"];
    method = new GTMethod::bayes_dirichlet(lambda); 
  }else{
    stop("Aggregation method not found.");
  }
  return method;
}

//[[Rcpp::export]]
List hclustcc_cpp(const List nb,const NumericMatrix& X,List method_obj) {
  
  
  
  
  // TODO collision detection a priori and priority queue as a map not multimap ? being consistent with hclust strategy for ties ?
  // TODO look at heller empirical bayes for prior specification
  
  
  int V = X.nrow();
  bool is_bayesian = method_obj.inherits("bayesian_gtmethod");
  int k_relaxed = 1;
  
  // compute data statistics needed for priors or distance
  GTMethod::GTMethod * method = init_method(method_obj);
  method->init(X);
  
  
  // data-structure creation
  // adjacency graph as an adjacency list
  std::vector<node> graph(2*V-1);
  // merge priority queue
  std::multimap<double,std::pair<int, int>,std::less<double>> priority_queue;
  // list of active nodes
  std::set<int> active_nodes;
  // current negative loglike
  double Llc = 0;
  double Pc  = 0;
  if(is_bayesian){
    Pc += -lgamma(V+1);
  }
  
  for(int i=0; i<nb.length(); ++i){
    if(nb[i]!=R_NilValue) {
      NumericVector nbi = as<NumericVector>(nb[i]);
      node cnode = node(method->init_node(i,X(i,_)));
      Llc+=cnode.height;
      active_nodes.insert(i);
      for(int n=0; n<nbi.length(); ++n){
        int j = nbi[n];
        if(i!=j){
          node vnode = method->init_node(j,X(j,_));
          double d = method->dist(&cnode,&vnode);
          cnode.neibs.insert(std::make_pair(j,d));
          if(i<j){
            priority_queue.insert(std::make_pair(d,std::make_pair(i,j)));
          }
          
        }
      }
      graph[i]=cnode;
    }
  }
  
  
  
  
  
  // Lets Merge !
  NumericMatrix merge(V-1,2);
  NumericVector height(V-1);
  NumericVector queue_size(V-1);
  NumericVector Ll(V);
  NumericVector Prior(V);
  Ll[0]=Llc;
  Prior[0]=Pc;
  for(int imerge=0;imerge<(V-1);imerge++){
    queue_size[imerge]=priority_queue.size();
    int K_before = V-imerge; 
    int node_id = V+imerge;
    auto best_merge = priority_queue.begin();
    
    
    // deal with isolated regions (no more merge possibles with contiguity constrains)
    if(best_merge==priority_queue.end()){
      
      
      // Relax contiguity constraints to build a complete hierarchy
      //priority_queue.clear(); // already empty
      k_relaxed = V-imerge;
      for(auto it = active_nodes.begin(); it != active_nodes.end(); ++it ) {
        int i = *it;
        graph[i].neibs.clear();
        for(auto it_nei = active_nodes.begin(); it_nei != active_nodes.end(); ++it_nei ){
          int j = *it_nei;
          if(i!=j){
            double d = method->dist(&graph[i],&graph[j]);
            graph[i].neibs.insert(std::make_pair(j,d));
            if(i<j){
              priority_queue.insert(std::make_pair(d,std::make_pair(i,j)));
            }
          }
        }
      }
      best_merge = priority_queue.begin();
    }
    
    std::pair<int,int> edge = best_merge->second; 
    int g = std::get<0>(edge);
    int h = std::get<1>(edge);
    node node_g = graph[g];
    node node_h = graph[h];
    
    height[imerge]=best_merge->first;
    Ll[imerge+1] = Llc+height[imerge];
    Prior[imerge+1] = Prior[imerge];
    if(is_bayesian){
      // a checker
      //Ll[imerge+1]+=log(K_before-1)-log(V-K_before+1);
      //Prior[imerge+1]+=log(K_before-1)-log(V-K_before+1);
      Prior[imerge+1]+=lgamma(node_h.size+node_g.size+1)-lgamma(node_h.size+1)-lgamma(node_g.size+1);
    }
    Llc = Ll[imerge+1];
    
    // strore merge move in hclust format with 1 based indices
    if(g<V){
      merge(imerge,0)=-(g+1);
    }else{
      merge(imerge,0)=g-V+1;
    }
    
    if(h<V){
      merge(imerge,1)=-(h+1);
    }else{
      merge(imerge,1)=h-V+1;
    }
    
    // update actives_nodes
    active_nodes.erase(active_nodes.find(g));
    active_nodes.erase(active_nodes.find(h));
    active_nodes.insert(node_id);
    
    
    // create a new node
    node new_node = node(method->merge(node_id,&node_g,&node_h,height[imerge]));
    
    
    // update the graph and priority queue
    for(auto nei_g = node_g.neibs.begin();nei_g!=node_g.neibs.end();nei_g++){
      
      int i = g;
      int j = nei_g->first;
      double v = nei_g->second;
      // old link deletion in priority_queue
      auto search = priority_queue.equal_range(v);
      for (auto s = search.first; s != search.second; ++s){
        std::pair<int,int> edge = s->second;
        if(std::get<0>(edge)==std::min(i,j) && std::get<1>(edge)==std::max(i,j)){
          priority_queue.erase(s);
          break;
        }
      }
      
      // old link deletion in graph
      graph[j].neibs.erase(i);
      // new links in graph
      if(j!=h){
        // distance calculation
        double d = method->dist(&new_node,&graph[j]);
        new_node.neibs.insert(std::make_pair(j,d));
        graph[j].neibs.insert(std::make_pair(node_id,d));
      }
      
      
    }
    
    
    for(auto nei_h = node_h.neibs.begin();nei_h!=node_h.neibs.end();nei_h++){
      
      int i = h;
      int j = nei_h->first;
      double v = nei_h->second;
      
      // old link deletion in priority_queue
      auto search = priority_queue.equal_range(v);
      for (auto s = search.first; s != search.second; ++s){
        std::pair<int,int> edge = s->second;
        if(std::get<0>(edge)==std::min(i,j) && std::get<1>(edge)==std::max(i,j)){
          priority_queue.erase(s);
          break;
        }
      }
      
      
      // old link deletion in graph
      graph[j].neibs.erase(i);
      
      // new links in graphs
      if(j!=g){
        double d = method->dist(&new_node,&graph[j]);
        new_node.neibs.insert(std::make_pair(j,d));
        graph[j].neibs.insert(std::make_pair(node_id,d));
      }
      
    }
    
    
    // add the newly created node
    graph[node_id]=new_node;
    
    // add the new possible merges in the priority queue
    for(auto nei = new_node.neibs.begin();nei!=new_node.neibs.end();nei++){
      double d = nei->second;
      double j = nei->first;
      priority_queue.insert(std::make_pair(d,std::make_pair(j,node_id)));
    }
    
  }
  // Export Centers
  NumericMatrix centers(V-1,X.ncol());
  for(int i=V;i<(2*V-1);i++){
    node cnode = graph[i];
    centers(i-V,_)=cnode.x;
  }
  CharacterVector ch = colnames(X);
  colnames(centers) = ch;
  delete method;
  List res = List::create(Named("merge",merge),
                          Named("Ll",Ll),
                          Named("Prior",Prior),
                          Named("queue_size",queue_size),
                          Named("data",X),
                          Named("centers",centers),
                          Named("k.relaxed",k_relaxed));
  return res;
}





//[[Rcpp::export]]
List bayesian_hclustcc_cpp(const List nb,const NumericMatrix& X,List method_obj) {
  
  
  // start cholmod
  cholmod_common c ;
  cholmod_start (&c) ; /* start CHOLMOD */
  
  
  int V = X.nrow();
  
  
  // compute data statistics needed for priors or distance
  GTMethod::GTMethod * method = init_method(method_obj);
  method->init(X);
  
  
  // data-structure creation
  // adjacency graph as an adjacency list
  std::vector<bayesian_node> graph(2*V-1);
  // merge priority queue
  std::multimap<double,std::pair<int, int>,std::less<double>> priority_queue;
  // list of active nodes
  std::set<int> active_nodes;
  // current negative loglike
  double Llc = 0;
  
  //Rcout << "init" << std::endl;
  int nblinks = 0;
  for(int i=0; i<nb.length(); ++i){
    if(nb[i]!=R_NilValue) {
      NumericVector nbi = as<NumericVector>(nb[i]);
      bayesian_node cnode(method->init_node(i,X(i,_)));
      cnode.lognbtree = 0;
      cnode.i_inter = i;
      Llc+=cnode.height;
      for(int n=0; n<nbi.length(); ++n){
        int j = nbi[n];
        if(i!=j){
          abstract_node vnode = method->init_node(j,X(j,_));
          double d = method->dist(&cnode,&vnode);
          multiedge e = multiedge(1,d);
          e.add_edge(std::make_pair(i,j));
          cnode.neibs.insert(std::make_pair(j,e));
          
          if(i<j){
            nblinks++;
            priority_queue.insert(std::make_pair(d,std::make_pair(i,j)));
          }
          
        }
      }
      graph[i]=cnode;
    }
  }
  
  // initialize inter cluster factorization
  cholmod_sparse * L = inter_to_sparse(graph,V,nblinks,&c);
  cholmod_factor * Linter = cholmod_analyze(L, &c) ; /* analyze */
  cholmod_factorize(L, Linter, &c) ; 
  double logdetfull = cholmod_tools_logdet(Linter);
  int * permutation = cholmod_tools_iPerm(((int *)Linter->Perm),V);
  int inter_pivot = permutation[V-1];
  
  // initialize intra clusters factorisation
  // we will use the same permuation for both.
  for(int r=0;r<V;r++){
    graph[r].intra_pivot = permutation[r];
    graph[r].intra_nodes.insert(graph[r].intra_nodes.begin(),permutation[r]);
    graph[r].i_inter = permutation[r];
    active_nodes.insert(permutation[r]);
  }
  cholmod_factor* Lintra = cholmod_allocate_factor(V,&c);
  
  // Lets Merge !
  NumericMatrix merge(V-1,2);
  NumericVector height(V-1);
  NumericVector queue_size(V-1);
  NumericVector Ll(V);
  NumericVector PriorIntra(V);
  NumericVector PriorInter(V);
  NumericVector PriorK(V);

  Ll[0]=Llc;
  PriorIntra[0]=0;
  //Rcout << "Main loop" << std::endl;
  for(int imerge=0;imerge<(V-1);imerge++){
    Rcout << "Merge N°" << imerge <<  std::endl;
    
    
    queue_size[imerge]=priority_queue.size();
    int K_before = V-imerge; 
    int node_id = V+imerge;
    auto best_merge = priority_queue.begin();
    
    
    
    std::pair<int,int> edge = best_merge->second; 
    int g = std::get<0>(edge);
    int h = std::get<1>(edge);
    

    
    bayesian_node node_g = graph[g];
    bayesian_node node_h = graph[h];
    
    
    
    // strore merge move in hclust format with 1 based indices
    if(g<V){
      merge(imerge,0)=-(g+1);
    }else{
      merge(imerge,0)=g-V+1;
    }
    
    if(h<V){
      merge(imerge,1)=-(h+1);
    }else{
      merge(imerge,1)=h-V+1;
    }
    

    // update actives_nodes
    active_nodes.erase(active_nodes.find(node_h.i_inter));

    // create a new node
    bayesian_node new_node(method->merge(node_id,&node_g,&node_h,height[imerge]));
    
    
    // update intra node graph
    // add node from g
    new_node.intra_nodes = node_g.intra_nodes;
    new_node.intra_nodes.splice(new_node.intra_nodes.end(),node_h.intra_nodes);
    
    // get the cutset g,h
    std::vector<std::pair<int, int>> cutset = node_g.neibs.at(h).edges;

    // build cholevky factorization of new node
    cholmod_tools_Lup_intra(Lintra,node_g.intra_pivot,node_h.intra_pivot,node_h.intra_pivot_edges,&cutset,permutation,&c);
    // update remaining pivot edges
    new_node.intra_pivot = node_g.intra_pivot;
    new_node.intra_pivot_edges =  cholmod_tools_pivotedgesup_intra(V,node_g.intra_pivot,node_g.intra_pivot_edges,&cutset,permutation);
    // store lognbtree
    new_node.lognbtree =  cholmod_tools_logdet_subset(Lintra,new_node.intra_nodes);

    Rcout << "lognbtree - intra :" << new_node.lognbtree << std::endl;
    PriorIntra[imerge+1] = PriorIntra[imerge];
    PriorIntra[imerge+1]+= new_node.lognbtree-node_g.lognbtree-node_h.lognbtree;
    
    

    
    
    
    
    height[imerge]=best_merge->first;
    Ll[imerge+1] = Llc-method->dist(&node_g,&node_h);
    Llc = Ll[imerge+1];
    
    // update the inter-graph and priority queue
    for(auto nei_g = node_g.neibs.begin();nei_g!=node_g.neibs.end();nei_g++){
      
      int i = g;
      int j = nei_g->first;
      
      
      
      multiedge uv = nei_g->second;
      // old link deletion in priority_queue
      double cheight = uv.height;
      
      auto search = priority_queue.equal_range(cheight);
      for (auto s = search.first; s != search.second; ++s){
        std::pair<int,int> edge = s->second;
        if(std::get<0>(edge)==std::min(i,j) && std::get<1>(edge)==std::max(i,j)){
          priority_queue.erase(s);
          break;
        }
      }
      
      // get the olds links
      multiedge old_links = graph[j].neibs.at(i);
      // old link deletion in graph
      graph[j].neibs.erase(i);
      // new links in graph
      if(j!=h){
        // distance calculation
        double d = method->dist(&new_node,&graph[j]);
        old_links.height=d;
        new_node.neibs.insert(std::make_pair(j,old_links));
        graph[j].neibs.insert(std::make_pair(node_id,old_links));
      }
      
      
    }
    
    
    for(auto nei_h = node_h.neibs.begin();nei_h!=node_h.neibs.end();nei_h++){
      
      int i = h;
      int j = nei_h->first;
      
      
      
      multiedge uv = nei_h->second;
      double cheight = uv.height;
      // old link deletion in priority_queue
      auto search = priority_queue.equal_range(cheight);
      for (auto s = search.first; s != search.second; ++s){
        std::pair<int,int> edge = s->second;
        if(std::get<0>(edge)==std::min(i,j) && std::get<1>(edge)==std::max(i,j)){
          priority_queue.erase(s);
          break;
        }
      }
      
      
      // old link deletion in graph
      multiedge old_links = graph[j].neibs.at(i);
      // old link deletion in graph
      graph[j].neibs.erase(i);
      
      // new links in graphs
      if(j!=g){
        // first case the links was not created by inspecting the neibs of g
        if(new_node.neibs.find(j)==new_node.neibs.end()){
          double d = method->dist(&new_node,&graph[j]);
          old_links.height=d;
          new_node.neibs.insert(std::make_pair(j,old_links));
          graph[j].neibs.insert(std::make_pair(node_id,old_links));
        }else{
          // the link was already created during the inpextion of g neibs
          new_node.neibs.at(j).merge_edges(old_links);
          graph[j].neibs.at(node_id).merge_edges(old_links);
        }
      }
      
    }
    
    
    
    // update i_inter
    new_node.i_inter = node_g.i_inter;
    // add the newly created node
    graph[node_id]=new_node;
    Rcout << "lognbtree - inter :";
    
    // compute Linter

    double lnbtreeinter = 0;
    if(node_h.i_inter!=inter_pivot){
      cholmod_rowdel(node_h.i_inter,NULL,Linter,&c);
      if(node_g.i_inter!=inter_pivot){
        std::map<int, int> links_ha = build_inter_links(&node_h,&graph);
        cholmod_sparse * links_toadd_inter = cholmod_tools_ltoadd_inter(V, node_g.i_inter, node_h.i_inter, inter_pivot,links_ha,&c);
        cholmod_updown(true,links_toadd_inter,Linter,&c);
        cholmod_free_sparse(&links_toadd_inter, &c);
        std::map<int, int> links_hd = build_inter_links(&node_h,&graph);
        cholmod_sparse * links_todel_inter = cholmod_tools_ltodel_inter(V, node_g.i_inter,node_h.i_inter, inter_pivot,links_hd,&c);
        cholmod_updown(false,links_todel_inter,Linter,&c);
        cholmod_free_sparse(&links_todel_inter, &c);
      }
      lnbtreeinter = cholmod_tools_logdet_intra(Linter,&active_nodes);
    }else{
      // g is the new pivot
      inter_pivot = node_g.i_inter;
      cholmod_rowdel(node_g.i_inter,NULL,Linter,&c);
      lnbtreeinter = cholmod_tools_logdet_intra(Linter,&active_nodes);
    }
    
    PriorInter[imerge+1] = lnbtreeinter;
    Rcout << lnbtreeinter  << std::endl;
    
    
    
    
    
    // add the new possible merges in the priority queue
    for(auto nei = new_node.neibs.begin();nei!=new_node.neibs.end();nei++){
      multiedge e = nei->second;
      double d = e.height;
      double j = nei->first;
      
      // update the merge cost to incorporate the prior
      // std::set<int> old_pivot_edges,const std::vector<std::pair<int, int>> * cutset,
      double cc = cut_cost(Lintra,&new_node,&(graph[j]),permutation,&c);
      
      double nd = d-cc;
      graph[node_id].neibs.at(j).height=nd;
      graph[j].neibs.at(node_id).height=nd;
      priority_queue.insert(std::make_pair(nd,std::make_pair(j,node_id)));
    }
    
    
    
    
    // to check
    PriorK[imerge+1]=PriorK[imerge]-log(K_before-1)+log(V-K_before+1)+log(K_before);
    
    
    
  }
  // intialize first value of prior inter // last value of prior intra // log nb of spanning tree
  PriorInter[0]=PriorIntra[V-1];
  
  // Export Centers
  NumericMatrix centers(V-1,X.ncol());
  for(int i=V;i<(2*V-1);i++){
    abstract_node cnode = graph[i];
    centers(i-V,_)=cnode.x;
  }
  CharacterVector ch = colnames(X);
  colnames(centers) = ch;
  delete method;
  List res = List::create(Named("merge",merge),
                          Named("Ll",Ll),
                          Named("PriorIntra",PriorIntra),
                          Named("PriorInter",PriorInter),
                          Named("PriorK",PriorK),
                          Named("queue_size",queue_size),
                          Named("data",X),
                          Named("centers",centers),
                          Named("k.relaxed",1));
  
  cholmod_finish (&c) ; /* finish CHOLMOD */
  
  free(permutation);
  return res;
}




