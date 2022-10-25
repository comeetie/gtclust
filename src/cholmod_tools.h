#ifndef CHOLMOD_TOOLS
#define CHOLMOD_TOOLS

#include "SuiteSparse_config/SuiteSparse_config.h"
#include "CHOLMOD/Include/cholmod.h"
// #include <Rcpp.h>
// using namespace Rcpp;


int * cholmod_tools_iPerm(int arr[], int size) {
  
  // to store element to index mappings
  int * iperm = (int*)malloc(sizeof(int) * size);
  for (int i = 0; i < size; i++){
    iperm[arr[i]] = i ;
  }
  return(iperm);
}

double cholmod_tools_logdet(cholmod_factor * L){
  double logdet = 0;
  for(int r=0;r<L->n;r++){
    
    int pd = ((int*)L->p)[r];
    double val = ((double*)L->x)[pd];
    if(val>0){
      logdet+= log(val); 
    }
  }
  return logdet;
}


double cholmod_tools_logdet_intra(cholmod_factor * L,const std::set<int> * subset){
  double logdet = 0;
  for(auto it=subset->begin();it!=subset->end();++it){
    int r = *it;
    int pd = ((int*)L->p)[r];
    double val = ((double*)L->x)[pd];
    if(val>0){
      logdet+= log(val); 
    }
  }
  return logdet;
}

template<template<typename, typename> typename Container, typename T, typename Allocator>
double cholmod_tools_logdet_subset(cholmod_factor * L,const Container<T, Allocator>& subset){
  double logdet = 0;
  for(const T& it : subset){
    int r = it;
    int pd = ((int*)L->p)[r];
    double val = ((double*)L->x)[pd];
    if(val>0){
      logdet+= log(val); 
    }
  }
  return logdet;
}


cholmod_sparse* cholmod_tools_linkstocol(int size,int pivot,int old_pivot,const std::vector<std::pair<int, int>> * cutset,   cholmod_common * com) {
  cholmod_sparse *S;
  
  int nb_links = cutset->size();
  S = cholmod_allocate_sparse(size,nb_links,2*nb_links,true,false,0,CHOLMOD_REAL, com);
  
  ((int*)S->p)[0] = 0;
  int u,v;
  int r = 0;
  for (auto it=cutset->begin();it<cutset->end();++it){
    if(it->first > it->second){
      u=it->second;
      v=it->first;
    }else{
      u=it->first;
      v=it->second;
    }
    ((int*)S->p)[r+1] =((int*)S->p)[r]+2;
    int nz=0;
    if(u!=pivot & u!=old_pivot){
      ((int*)S->i)[r*2] = u;
      ((double*)S->x)[r*2] = 1;
      nz++;
    }
    if(v!=pivot & v!=old_pivot){
      ((int*)S->i)[r*2+nz] = v;
      ((double*)S->x)[r*2+nz] = -1;
      nz++;
    }
    ((int*)S->nz)[r]=nz;
    r++;
  }
  
  return(S);
}


// we pass intra_old_pivot_edges by value, to avoid any side effect
cholmod_sparse * cholmod_tools_oldpivotcol(int size,int pivot, int old_pivot,std::set<int> intra_oldpivot_edges,const std::vector<std::pair<int, int>>  * cutset, cholmod_common * com){
  // find any link in the cutset that deal with old_pivot
  for (auto it=cutset->begin();it<cutset->end();++it){
    if(it->first==old_pivot){
      intra_oldpivot_edges.insert(it->second);
    }
    if(it->second==old_pivot){
      intra_oldpivot_edges.insert(it->first);
    }
  }

  // number of links
  int nb_links = intra_oldpivot_edges.size();
  // if their is a link with the pivot erase it
  auto it = intra_oldpivot_edges.find(pivot);
  if(it!=intra_oldpivot_edges.end()){
    intra_oldpivot_edges.erase(it);
  }
  // insert the diagonal 
  intra_oldpivot_edges.insert(old_pivot);
  
  // create the cholmod sparse column
  cholmod_sparse * new_col = cholmod_allocate_sparse(size,1,intra_oldpivot_edges.size(),true,true,0,CHOLMOD_REAL, com);
  int r=0;
  for (auto it=intra_oldpivot_edges.begin();it!=intra_oldpivot_edges.end();++it){
    int v = *it;
    ((int*)new_col->i)[r] = v;
    if(v==old_pivot){
      ((double*)new_col->x)[r] = nb_links;
    }else{
      ((double*)new_col->x)[r] = -1;
    }
    r++;
  }
  ((int*)new_col->p)[0]=0;
  ((int*)new_col->p)[1]=r;
  return(new_col);
  
}




cholmod_factor * cholmod_tools_Lup_intra(cholmod_factor * F,int pivot, int old_pivot,std::set<int> old_pivot_edges,const std::vector<std::pair<int, int>> * cutset,int * permutation, cholmod_common * com){
  
  std::vector<std::pair<int, int>> cutset_loc(cutset->size());
  int i=0;
  for (auto it=cutset->begin(); it !=cutset->end(); ++it){
    cutset_loc[i] = std::make_pair(permutation[it->first],permutation[it->second]);
    i++;
  }

  cholmod_sparse* links1 = cholmod_tools_linkstocol(F->n,pivot,old_pivot,&cutset_loc,com);

  cholmod_updown(true,links1,F,com);
  
  cholmod_free_sparse (&links1, com);

  cholmod_sparse * col1 = cholmod_tools_oldpivotcol(F->n,pivot,old_pivot,old_pivot_edges,&cutset_loc,com);

  cholmod_rowadd(old_pivot,col1,F,com);
  
  cholmod_free_sparse (&col1, com);
  
  return (F) ;
}

double cholmod_tools_cutcost(cholmod_factor * F,const bayesian_node * node_g,const bayesian_node * node_h,std::set<int> old_pivot_edges,const std::vector<std::pair<int, int>> * cutset,int * permutation, cholmod_common * com){
  
  int pivot = node_g->intra_pivot; 
  int old_pivot=node_h->intra_pivot;
  std::vector<std::pair<int, int>> cutset_loc(cutset->size());
  int i=0;
  for (auto it=cutset->begin(); it !=cutset->end(); ++it){
    cutset_loc[i] = std::make_pair(permutation[it->first],permutation[it->second]);
    i++;
  }
  double lognbtree=0;
  // cholmod_sparse* links1 = cholmod_tools_linkstocol(F->n,pivot,old_pivot,&cutset_loc,com);
  // cholmod_updown(true,links1,F,com);
  // cholmod_sparse * col1 = cholmod_tools_oldpivotcol(F->n,pivot,old_pivot,old_pivot_edges,&cutset_loc,com);
  // cholmod_rowadd(old_pivot,col1,F,com);
  // 
  // // computelogdet
  // double lognbtree =  cholmod_tools_logdet_subset(F,node_g->intra_nodes)+cholmod_tools_logdet_subset(F,node_h->intra_nodes);
  // 
  // // revert the changes
  // cholmod_rowdel(old_pivot,NULL,F,com);
  // cholmod_free_sparse (&col1, com);
  // cholmod_updown(false,links1,F,com);
  // cholmod_free_sparse (&links1, com);
  
  return (lognbtree) ;
}




std::set<int> cholmod_tools_pivotedgesup_intra(int size,int pivot,std::set<int>  intra_pivot_edges, const std::vector<std::pair<int, int>> * cutset,int * permutation){
  for (auto it=cutset->begin();it<cutset->end();++it){
    if(permutation[it->first]==pivot){
      intra_pivot_edges.insert(permutation[it->second]);
    }
    if(permutation[it->second]==pivot){
      intra_pivot_edges.insert(permutation[it->first]);
    }
  }
  return intra_pivot_edges;
}




cholmod_sparse* cholmod_tools_colinter(std::map<int, int> inter,int i, int n, int pivot,cholmod_common* com){
  auto search = inter.find(pivot);
  int nbval = inter.size();
  if(search!=inter.end()){
    inter.erase(search);
    nbval--;
  }
  cholmod_sparse * col = cholmod_allocate_sparse(n,1,nbval,true,true,0,CHOLMOD_REAL, com);
  int vp = 0;
  for (auto itr = inter.begin(); itr!=inter.end();++itr){
    ((int *)col->i)[vp]=itr->first;
    ((double *)col->x)[vp]=itr->second;
    vp++;
  }
  ((int *)col->p)[0]=0;
  ((int *)col->p)[1]=nbval;
  return col;
}




cholmod_sparse* cholmod_tools_ltoadd_inter(int n, int g, int h, int pivot, std::map<int, int> inter_h,cholmod_common* com ) {
  
  auto search = inter_h.find(g);
  //remove g-h from h 
  if(search!=inter_h.end()){
    inter_h.erase(search);
  }
  //remove diagonals
  search = inter_h.find(h);
  if(search!=inter_h.end()){
    inter_h.erase(search);
  }
  // links with the pivot ?
  search = inter_h.find(pivot);
  double piv_val = 0;
  int link_piv = 0;
  if(search!=inter_h.end()){
    link_piv = 1;
    piv_val = - search->second;
    inter_h.erase(search);
  }
  
  int nb_links = inter_h.size();
  
  cholmod_sparse *S;
  
  S = cholmod_allocate_sparse(n,nb_links+link_piv,2*nb_links+link_piv,true,true,0,CHOLMOD_REAL, com);
  int r = 0;
  int u,v;
  ((int*)S->p)[0] = 0;
  // ! handle case were h or h is the pivot not done
  for (auto itr = inter_h.begin(); itr!=inter_h.end();++itr){
    ((int*)S->p)[r+1] =((int*)S->p)[r]+2;
    if(itr->first>g){
      u=g;
      v=itr->first;
    }else{
      u=itr->first;
      v=g;
    }
    ((int*)S->i)[r*2] = u;
    ((double*)S->x)[r*2] = -sqrt(-itr->second);
    ((int*)S->i)[r*2+1] = v;
    ((double*)S->x)[r*2+1] = sqrt(-itr->second);
    r++;
  }
  if(link_piv==1){
    ((int*)S->p)[r+1] =((int*)S->p)[r]+1;
    ((int*)S->i)[r*2] = g;
    ((double*)S->x)[r*2] = sqrt(piv_val);
  }
  
  //print_sparse(S);
  return(S);
}



cholmod_sparse* cholmod_tools_ltodel_inter(int n, int g, int h, int pivot, std::map<int, int> inter_h,cholmod_common* com ) {
  
  //remove diagonals
  auto search = inter_h.find(h);
  if(search!=inter_h.end()){
    inter_h.erase(search);
  }
  
  // links with pivot ?
  search = inter_h.find(pivot);
  if(search!=inter_h.end()){
    inter_h.erase(search);
  }
  
  
  
  int nb_links = inter_h.size();
  cholmod_sparse *S;
  
  S = cholmod_allocate_sparse(n,nb_links,nb_links,true,true,0,CHOLMOD_REAL, com);
  int r = 0;
  ((int*)S->p)[0] = 0;
  for (auto itr = inter_h.begin(); itr!=inter_h.end();++itr){
    ((int*)S->p)[r+1] =((int*)S->p)[r]+1;
    ((int*)S->i)[r] = itr->first;
    ((double*)S->x)[r] = sqrt(-itr->second);
    r++;
  }
  //print_sparse(S);
  return(S);
}


cholmod_sparse * inter_to_sparse(std::vector<bayesian_node> graph,int V,int nblinks,cholmod_common * com){
  
  cholmod_sparse * L = cholmod_allocate_sparse(V,V,nblinks+V,true,false,-1,CHOLMOD_REAL, com);
  int ichol = 0;
  
  ((int*)L->p)[0]=0;
  for(int icol=0;icol<V;icol++) {
    std::map<int, int> wlinks;
    int diagval = 0;
    for (auto itr = graph[icol].neibs.begin(); itr!=graph[icol].neibs.end();++itr){
      int jg = itr->first;
      int j = graph[jg].i_inter;
      multiedge e = itr->second;
      diagval+=e.size;
      if(j>icol){
        wlinks.insert(std::make_pair(j,-e.size));
      }
    }
    wlinks.insert(std::make_pair(icol,diagval));
    int nbl = wlinks.size();
    ((int*)L->p)[icol+1]=((int*)L->p)[icol]+nbl;
    if(icol!=(V-1)){
      for(auto itl =wlinks.begin();itl!=wlinks.end();++itl){
        //Rcout << icol << "--" << itl-> first << itl->second << std::endl;
        if(itl->first!=(V-1)){
          ((int*)L->i)[ichol]=itl->first;
          ((double*)L->x)[ichol]=itl->second;
        }else{
          nbl--;
        }
        ichol++;
      }
      ((int*)L->nz)[icol]=nbl;
    }else{
      ((int*)L->i)[ichol]=icol;
      ((double*)L->x)[ichol]=1;
      ((int*)L->nz)[icol]=1;
      ((int*)L->p)[icol+1]=((int*)L->p)[icol]+1;
      ichol++;
    }
    
  }
  //print_sparse(L);
  return(L);
}


// void print_factor(cholmod_factor * L){
//   Rcout << "is ll :" << L->is_ll << std::endl;
//   Rcout << "is monotonic :" << L->is_monotonic << std::endl;
//   Rcout << "is super :" << L->is_super << std::endl;
//   Rcout << "IX: " ;
//   for(int v=0;v<L->nzmax;v++){
//     Rcout << "(" << ((int*)L->i)[v] << ", " << ((double*)L->x)[v] << ")";
//   }
//   Rcout << std::endl;
//   Rcout << "PNZ: ";
//   for(int v=0;v<L->n;v++){
//     Rcout << "(" << ((int*)L->p)[v] << ", " << ((int*)L->nz)[v] << ")";
//   }
//   Rcout << std::endl;
//   Rcout << "NP: ";
//   for(int v=0;v<(L->n+2);v++){
//     Rcout << "(" << ((int*)L->prev)[v] << ", " << ((int*)L->next)[v] << ")";
//   }
//   Rcout << std::endl;
// }
// 
// 
// void print_sparse(cholmod_sparse * S){
//   Rcout << "is packed :" << S->packed << std::endl;
//   Rcout << "IX: " ;
//   for(int v=0;v<S->nzmax;v++){
//     Rcout << "(" << ((int*)S->i)[v] << ", " << ((double*)S->x)[v] << ")";
//   }
//   Rcout << std::endl;
//   Rcout << "P: ";
//   for(int v=0;v<(S->ncol+1);v++){
//     Rcout << ((int*)S->p)[v] << ", ";
//   }
//   Rcout << std::endl;
//   if(!S->packed){
//     Rcout << "NZ: ";
//     for(int v=0;v<(S->ncol);v++){
//       Rcout << ((int*)S->nz)[v] << ", ";
//     }
//     Rcout << std::endl;
//   }
// }


#endif

