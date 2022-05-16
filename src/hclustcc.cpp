#include <Rcpp.h>
#include "GTMethod.h"
#include "node.h"
using namespace Rcpp;




void print_pq(std::multimap<float,std::pair<int, int>,std::less<float>> priority_queue){
  for(auto it=priority_queue.begin();it!=priority_queue.end();it++){
    std::pair<int,int> edge = it->second;
    Rcout << edge.first << " - " << edge.second << ":" << it->first << std::endl;
  }
}

void print_pqhead(std::multimap<float,std::pair<int, int>,std::less<float>> priority_queue){
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
    float beta = method_obj["beta"];
    method = new GTMethod::bayes_mom(beta); 
  }
  return method;
}

//[[Rcpp::export]]
List hclustcc_cpp(const List nb,const NumericMatrix X,List method_obj) {

  
  
  // TODO specific class for bayesian results optimal k* and test value results, height definition ? 
  // TODO collision detection a priori and priority queue as a map not multimap ? being consistent with hclust strategy for ties ?
  // TODO isolated island completion 
  // TODO interface // prior specification 
  // TODO bayesdgmm method ?
  // TODO look at heller empirical bayes for prior specification
  // TODO test and check chi2 
  
  int V = X.nrow();

  // cst added in case of collision
  // float collision_eps=1e-12;
  
  // compute data statistics needed for priors or distance
  GTMethod::GTMethod * method = init_method(method_obj);
  method->init(X);

  
  // data-structure creation
  std::vector<node> graph(2*V-1);
  std::multimap<float,std::pair<int, int>,std::less<float>> priority_queue;
  for(int i=0; i<nb.length(); ++i){
    if(nb[i]!=R_NilValue) {
      NumericVector nbi = as<NumericVector>(nb[i]);
      node cnode = method->init_node(i,X(i,_));
      for(int n=0; n<nbi.length(); ++n){
        int j = nbi[n];
        if(i!=j){
          node vnode = method->init_node(j,X(j,_));
          float d = method->dist(cnode,vnode);
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
  
  for(int imerge=0;imerge<(V-1);imerge++){
    

    int node_id = V+imerge;
    auto best_merge = priority_queue.begin();
    
    
    // deal with isolated regions (no more merge possible and node)
    if(best_merge==priority_queue.end()){

      // Rcout << "isolated regions" << std::endl;
      for(int nb_miss_merge=node_id;nb_miss_merge<(V-1);nb_miss_merge++){
        merge(nb_miss_merge,1)=0;
        merge(nb_miss_merge,2)=0;
      }
      break;
    }
    
    std::pair<int,int> edge = best_merge->second; 
    int g = std::get<0>(edge);
    int h = std::get<1>(edge);
    node node_g = graph[g];
    node node_h = graph[h];
    
    height[imerge]=best_merge->first;
    
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
    
    // create a new node
    node new_node = method->merge(node_id,node_g,node_h,height[imerge]);
    
    // update the graph and priority queue
    for(auto nei_g = node_g.neibs.begin();nei_g!=node_g.neibs.end();nei_g++){
      
      int i = g;
      int j = nei_g->first;
      float v = nei_g->second;
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
        float d = method->dist(new_node,graph[j]);
        new_node.neibs.insert(std::make_pair(j,d));
        graph[j].neibs.insert(std::make_pair(node_id,d));
      }
      
      
    }
    
    
    for(auto nei_h = node_h.neibs.begin();nei_h!=node_h.neibs.end();nei_h++){
      
      int i = h;
      int j = nei_h->first;
      float v = nei_h->second;

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
        
        float d = method->dist(new_node,graph[j]);
        new_node.neibs.insert(std::make_pair(j,d));
        graph[j].neibs.insert(std::make_pair(node_id,d));
      }
      
    }
    
    
    // add the newly created node
    graph[node_id]=new_node;
    
    // add the new possible merges in the priority queue
    for(auto nei = new_node.neibs.begin();nei!=new_node.neibs.end();nei++){
      float d = nei->second;
      float j = nei->first;
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

  return List::create(Named("merge",merge),Named("height",height),Named("data",X),Named("centers",centers));
  
}


