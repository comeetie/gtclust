#include <Rcpp.h>
#include <RcppEigen.h>
extern "C" {
#include <RcppEigenCholmod.h>
}
#include "GTMethod.h"
#include "node.h"
#include "LaplaceDirichlet.h"
using namespace Rcpp;


//[[Rcpp::export]]
Eigen::SparseMatrix<double> sparseBdiag(Rcpp::List B_list)
{
  int K = B_list.length();
  Eigen::VectorXi B_cols(K);
  for(int k=0;k<K;k++) {
    Eigen::SparseMatrix<double> Bk = B_list(k);
    B_cols[k] = Bk.cols();
  }
  int sumCols = B_cols.sum();
  Eigen::SparseMatrix<double> A(sumCols,sumCols);
  int startCol=0, stopCol=0, Bk_cols;
  for(int k=0;k<K;k++) {
    Eigen::SparseMatrix<double> Bk = B_list(k);
    Bk_cols = Bk.cols();
    stopCol = startCol + Bk_cols;
    for(int j=startCol;j<stopCol;j++){
      A.startVec(j);
      for(Eigen::SparseMatrix<double,0,int>::InnerIterator it(Bk,j-startCol); it; ++it) {
        A.insertBack(it.row()+startCol,j) = it.value();
      }
    }
    startCol = stopCol;
  }
  A.finalize();
  return A;
}


//[[Rcpp::export]]
Eigen::SparseMatrix<double> remove1row1col(Eigen::SparseMatrix<double> L){
  int nbc = L.cols();
  Eigen::SparseMatrix<double> Lr(nbc-1,nbc-1);
  Eigen::VectorXi csizes(nbc-1);
  for(int j=0;j<(nbc-1);j++){
    csizes(j)=L.col(j+1).nonZeros();
  }
  for(Eigen::SparseMatrix<double>::InnerIterator it(L,0); it; ++it) {
    if(it.row()>0){
      csizes(it.row()-1)+=-1;
    }
  }
  Lr.reserve(csizes);
  for(int j=1;j<nbc;j++){
    for(Eigen::SparseMatrix<double>::InnerIterator it(L,j); it; ++it) {
      if(it.row()>0){
        Lr.insert(it.row()-1,j-1) = it.value();
      }
    }
  }
  Lr.makeCompressed();
  return Lr;
}
  
//[[Rcpp::export]]
Eigen::SparseMatrix<double> buildLaplacianR(Eigen::SparseMatrix<double> Lg, Eigen::SparseMatrix<double> Lh, NumericMatrix cutset,NumericVector permutation)
{
  int nbc_g = Lg.cols();
  int nbc_h = Lh.cols();
  Eigen::SparseMatrix<double> L(nbc_g+nbc_h,nbc_g+nbc_h);
  Eigen::VectorXi csizes(nbc_g+nbc_h);
  for(int j=0;j<nbc_g;j++){
    csizes(j)=Lg.col(j).nonZeros();
  }
  for(int j=0;j<nbc_h;j++){
    csizes(j+nbc_g)=Lh.col(j).nonZeros();
  }
  for (int ir=0; ir <cutset.cols(); ++ir){
    int i = permutation[cutset(ir,0)];
    int j = permutation[cutset(ir,1)];
    csizes(i) += 1;
    csizes(j) += 1;
  }

  L.reserve(csizes);
    
    
  
  // insert Lg
  for(int j=0;j<nbc_g;j++){
    for(Eigen::SparseMatrix<double,0,int>::InnerIterator it(Lg,j); it; ++it) {
      L.insert(it.row(),j) = it.value();
    }
  }
  
  // insert Lh
  for(int j=nbc_g;j<(nbc_g+nbc_h);j++){
    for(Eigen::SparseMatrix<double,0,int>::InnerIterator it(Lh,j-nbc_g); it; ++it) {
      L.insert(it.row()+nbc_g,j) = it.value();
    }
  }
  
  // insert the cutset
  for (int ir=0; ir <cutset.cols(); ++ir){
    int i = permutation[cutset(ir,0)];
    int j = permutation[cutset(ir,1)];
    L.coeffRef(i,i)+=1;
    L.coeffRef(j,j)+=1;
    L.insert(i,j) = -1;
    L.insert(j,i) = -1;
  }
  //L.finalize();
  L.makeCompressed();
  return L;
}

Eigen::SparseMatrix<double> buildLaplacian(Eigen::SparseMatrix<double> Lg, Eigen::SparseMatrix<double> Lh, std::vector<std::pair<int, int>> cutset,NumericVector permutation)
{
  int nbc_g = Lg.cols();
  int nbc_h = Lh.cols();
  Eigen::SparseMatrix<double> L(nbc_g+nbc_h,nbc_g+nbc_h);
  Eigen::VectorXi csizes(nbc_g+nbc_h);
  for(int j=0;j<nbc_g;j++){
    csizes(j)=Lg.col(j).nonZeros();
  }
  for(int j=0;j<nbc_h;j++){
    csizes(j+nbc_g)=Lh.col(j).nonZeros();
  }
  for (auto it=cutset.begin(); it !=cutset.end(); ++it){
    int i = permutation[it->first];
    int j = permutation[it->second];
    csizes(i) += 1;
    csizes(j) += 1;
  }
  
  L.reserve(csizes);
  
  
  // insert Lg
  for(int j=0;j<nbc_g;j++){
    for(Eigen::SparseMatrix<double,0,int>::InnerIterator it(Lg,j); it; ++it) {
      L.insert(it.row(),j) = it.value();
    }
  }
  
  // insert Lh
  for(int j=nbc_g;j<(nbc_g+nbc_h);j++){
    for(Eigen::SparseMatrix<double,0,int>::InnerIterator it(Lh,j-nbc_g); it; ++it) {
      L.insert(it.row()+nbc_g,j) = it.value();
    }
  }
  
  // insert the cutset
  for (auto it=cutset.begin(); it !=cutset.end(); ++it){
    int i = permutation[it->first];
    int j = permutation[it->second];
    L.coeffRef(i,i)+=1;
    L.coeffRef(j,j)+=1;
    L.insert(i,j) = -1;
    L.insert(j,i) = -1;
  }
  return L;
}

void print_pq(std::multimap<double,std::pair<int, int>,std::less<double>> priority_queue){
  for(auto it=priority_queue.begin();it!=priority_queue.end();it++){
    std::pair<int,int> edge = it->second;
    Rcout << edge.first << " - " << edge.second << ":" << it->first << std::endl;
  }
}

void print_pqhead(std::multimap<double,std::pair<int, int>,std::less<double>> priority_queue){
  int i =0;
  for(auto it=priority_queue.begin();it!=priority_queue.end();it++){
    std::pair<int,int> edge = it->second;
    Rcout << edge.first << " - " << edge.second << ":" << it->first << std::endl;
    if(i>10){
      break;
    }
    i++;
  }
}


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
  
  // Eigen solver
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver;



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
  
  
  for(int i=0; i<nb.length(); ++i){
    if(nb[i]!=R_NilValue) {
      NumericVector nbi = as<NumericVector>(nb[i]);
      bayesian_node cnode(method->init_node(i,X(i,_)));
      cnode.intra_nodes.push_back(i);
      Eigen::SparseMatrix<double> Lc(1,1);
      cnode.Laplacian = Lc;
      cnode.lognbtree = 0;
      Llc+=cnode.height;
      active_nodes.insert(i);
      for(int n=0; n<nbi.length(); ++n){
        int j = nbi[n];
        if(i!=j){
          abstract_node vnode = method->init_node(j,X(j,_));
          double d = method->dist(&cnode,&vnode);
          multiedge e = multiedge(1,d);
          e.add_edge(std::make_pair(i,j));
          cnode.neibs.insert(std::make_pair(j,e));
          
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
  NumericVector PriorIntra(V);
  NumericVector permutation(V);
  Ll[0]=Llc;
  PriorIntra[0]=0;
  for(int imerge=0;imerge<(V-1);imerge++){
    queue_size[imerge]=priority_queue.size();
    int K_before = V-imerge; 
    int node_id = V+imerge;
    auto best_merge = priority_queue.begin();
    
    
    
    std::pair<int,int> edge = best_merge->second; 
    int g = std::get<0>(edge);
    int h = std::get<1>(edge);
    
    
    
    bayesian_node node_g = graph[g];
    bayesian_node node_h = graph[h];
    
    
    
    
    height[imerge]=best_merge->first;
    Ll[imerge+1] = Llc+height[imerge];
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
    bayesian_node new_node(method->merge(node_id,&node_g,&node_h,height[imerge]));
    
    // update intra node graph
    // add node from g
    new_node.intra_nodes.insert(new_node.intra_nodes.end(),node_g.intra_nodes.begin(),node_g.intra_nodes.end());
    // add nodes from h and change permutation vector for them
    int offset = node_g.size;
    for (auto it=node_h.intra_nodes.begin(); it !=node_h.intra_nodes.end(); ++it){
      new_node.intra_nodes.push_back(*it);
      permutation[*it]=permutation[*it]+offset;
    }
    // get the cutset g,h
    std::vector<std::pair<int, int>> cutset = node_g.neibs.at(h).edges;
    
    // build le laplacian
    new_node.Laplacian =  buildLaplacian(node_g.Laplacian,node_h.Laplacian,cutset,permutation);
    // if(new_node.Laplacian.cols()<15){
    //   Rcout << Eigen::MatrixXd(new_node.Laplacian) << std::endl;
    // }

    // // Eigen solver
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    // Compute the numerical factorization 
    solver.compute(remove1row1col(new_node.Laplacian)); 
    cholmod_sparse Lchol = Eigen::viewAsCholmod(remove1row1col(new_node.Laplacian));
    


    if(solver.info()!=Eigen::Success) {
      stop("Eigen decomp failed");
    }else{
      double ldet = solver.logAbsDeterminant();
      new_node.lognbtree=0;
      if(ldet>1e-10){
        new_node.lognbtree= ldet;
      }
      
      PriorIntra[imerge+1] = PriorIntra[imerge];
      PriorIntra[imerge+1]+= new_node.lognbtree-node_g.lognbtree-node_h.lognbtree;
    }

    
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
    
    
    // add the newly created node
    graph[node_id]=new_node;
    
    // add the new possible merges in the priority queue
    for(auto nei = new_node.neibs.begin();nei!=new_node.neibs.end();nei++){
      multiedge e = nei->second;
      double d = e.height;
      double j = nei->first;
      priority_queue.insert(std::make_pair(d,std::make_pair(j,node_id)));
    }
    
  }
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
                          Named("queue_size",queue_size),
                          Named("data",X),
                          Named("centers",centers),
                          Named("permutation",permutation));
  return res;
}






