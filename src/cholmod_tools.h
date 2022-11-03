#ifndef CHOLMOD_TOOLS
#define CHOLMOD_TOOLS

#include "SuiteSparse_config/SuiteSparse_config.h"
#include "CHOLMOD/Include/cholmod.h"
#include "node.h"
#include <Rcpp.h>
using namespace Rcpp;

void print_factor(cholmod_factor * L){
  Rcout << "is ll :" << L->is_ll << std::endl;
  Rcout << "is monotonic :" << L->is_monotonic << std::endl;
  Rcout << "is super :" << L->is_super << std::endl;
  Rcout << "IX: " ;
  for(int v=0;v<L->nzmax;v++){
    Rcout << "(" << ((int*)L->i)[v] << ", " << ((double*)L->x)[v] << ")";
  }
  Rcout << std::endl;
  Rcout << "PNZ: ";
  for(int v=0;v<L->n;v++){
    Rcout << "(" << ((int*)L->p)[v] << ", " << ((int*)L->nz)[v] << ")";
  }
  Rcout << std::endl;
  Rcout << "NP: ";
  for(int v=0;v<(L->n+2);v++){
    Rcout << "(" << ((int*)L->prev)[v] << ", " << ((int*)L->next)[v] << ")";
  }
  Rcout << std::endl;
}

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
    }else{
      stop("Error: Cholevshy factorization update failed");
    }
  }
  return logdet;
}


double cholmod_tools_logdet_subset_delta(cholmod_factor * L,const std::set<int,std::greater<int>>& subset,double Lintra_diag[],int imin){
  double logdet_delta = 0;
  for(auto it=subset.begin();it!=subset.end();++it){
    int r = *it;
    if(r>=imin){
      int pd = ((int*)L->p)[r];
      double val = ((double*)L->x)[pd];
      if(val>0){
        logdet_delta+= log(val)-log(Lintra_diag[r]);
        Lintra_diag[r]=val;
      }else{
        logdet_delta-=log(Lintra_diag[r]);
        Lintra_diag[r]=1;
        Rcout << "warning:" << r << " : " << val << std::endl;
        //stop("Error: Cholevshy factorization update failed");
      }
    }else{
      break;
    }
  }
  return logdet_delta;
}


double cholmod_tools_logdet_subset_delta_nochange(cholmod_factor * L,const std::set<int,std::greater<int>>& subset,double Lintra_diag[],int imin){
  double logdet_delta = 0;
  for(auto it=subset.begin();it!=subset.end();++it){
    int r = *it;
    if(r>=imin){
      int pd = ((int*)L->p)[r];
      double val = ((double*)L->x)[pd];
      if(val>0){
        logdet_delta+= log(val)-log(Lintra_diag[r]);
      }else{
        logdet_delta-= log(Lintra_diag[r]);
        Rcout << "warning:" << r << " : " << val << std::endl;
        //stop("Error: Cholevshy factorization update failed");
      }
    }else{
      break;
    }
  }
  return logdet_delta;
}

double cholmod_tools_logdet_subset(cholmod_factor * L,const std::set<int,std::greater<int>>& subset,int imin){
  double logdet = 0;
  for(auto it=subset.begin();it!=subset.end();++it){
    int r = *it;
    if(r>=imin){
      int pd = ((int*)L->p)[r];
      double val = ((double*)L->x)[pd];
      if(val>0){
        logdet+= log(val); 
      }else{
        stop("Error: Cholevshy factorization update failed");
        //Rcout << "warning:" << r << " : " << val << std::endl;
      }
    }else{
      break;
    }
  }
  return logdet;
}

double cholmod_tools_logdet_subset(cholmod_factor * L,const std::set<int,std::greater<int>>& subset){
  double logdet = 0;
  for(auto it=subset.begin();it!=subset.end();++it){
    int r = *it;
    int pd = ((int*)L->p)[r];
    double val = ((double*)L->x)[pd];
    if(val>0){
      logdet+= log(val); 
    }else{
      stop("Error: Cholevshy factorization update failed");
      //Rcout << "warning:" << r << " : " << val << std::endl;
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
    if(u!=pivot && u!=old_pivot){
      ((int*)S->i)[r*2] = u;
      ((double*)S->x)[r*2] = 1;
      nz++;
    }
    if(v!=pivot && v!=old_pivot){
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




void cholmod_tools_Lup_intra(cholmod_factor * F,int pivot, int old_pivot,std::set<int> old_pivot_edges,const std::vector<std::pair<int, int>> * cutset,int * permutation, cholmod_common * com){
  
  std::vector<std::pair<int, int>> cutset_loc(cutset->size());
  int i=0;
  int imin = old_pivot;
  for (auto it=cutset->begin(); it !=cutset->end(); ++it){
    if(permutation[it->first]<imin){
      imin=permutation[it->first];
    }
    if(permutation[it->second]<imin){
      imin=permutation[it->second];
    }
    cutset_loc[i] = std::make_pair(permutation[it->first],permutation[it->second]);
    i++;
  }
  
  cholmod_sparse* links1 = cholmod_tools_linkstocol(F->n,pivot,old_pivot,&cutset_loc,com);

  cholmod_updown(true,links1,F,com);
  

  
  cholmod_free_sparse (&links1, com);

  
  for (auto it=old_pivot_edges.begin();it!=old_pivot_edges.end();++it){
    if(*it<imin){
      imin = *it; 
    }
  }

  cholmod_sparse * col1 = cholmod_tools_oldpivotcol(F->n,pivot,old_pivot,old_pivot_edges,&cutset_loc,com);

  cholmod_rowadd(old_pivot,col1,F,com);
  

  cholmod_free_sparse (&col1, com);
  
}


subfactor * cholmod_tools_save_subfactor(cholmod_factor * F,const std::set<int,std::greater<int>> & node_g_intra_nodes,int iming,const std::set<int,std::greater<int>> & node_h_intra_nodes, int iminh,cholmod_common * com){
  
  // get the sizes
  int tot_nz=0;
  int tot_col=0;
  for(auto it=node_g_intra_nodes.begin();it!=node_g_intra_nodes.end();++it){
    int ic =*it;
    if(ic>=iming){
      tot_nz += ((int*)F->nz)[ic];
      tot_col++;
    }else{
      break;
    }
  }
  for(auto it=node_h_intra_nodes.begin();it!=node_h_intra_nodes.end();++it){
    int ic =*it;
    if(ic>=iminh){
      tot_nz += ((int*)F->nz)[ic];
      tot_col++;
    }else{
      break;
    }
  }
  //Rcout << "tot_col :" << tot_col << std::endl;
  //Rcout << "tot_nz :" << tot_nz << std::endl;
  int  * il = (int*)malloc((tot_nz+1)*sizeof(int));
  double  * xl = (double*)malloc((tot_nz+1)*sizeof(double));
  int  * p = (int*)malloc((tot_col+1)*sizeof(int));
  int * cols = (int*)malloc((tot_col+1)*sizeof(int));
  int * nbz = (int*)malloc((tot_col+1)*sizeof(int));
  
  // int  il[50000];
  // double  xl[50000];
  // int  p[5000];
  // int cols[5000];
  // int nbz[5000];
  int ival=0;
  int ccol=0;
  for(auto it=node_g_intra_nodes.begin();it!=node_g_intra_nodes.end();++it){
    int ic = *it;
    cols[ccol] = ic;
    p[ccol]= ival;
    
    if(ic>=iming){
      // copy all col values
      int pd = ((int*)F->p)[ic];
      int cnz = ((int*)F->nz)[ic];
      nbz[ccol]= cnz;
      for (int i=pd;i<pd+cnz;++i){
        // Rcout << ival << " : (" << ((int*)F->i)[i] << "," << ((double*)F->x)[i] << std::endl;
        il[ival] = ((int*)F->i)[i];
        xl[ival] = ((double*)F->x)[i];
        ival++;
      }
      ccol++;
    }else{
      break;
    }
  }
  
  for(auto it=node_h_intra_nodes.begin();it!=node_h_intra_nodes.end();++it){
    int ic = *it;
    cols[ccol] = ic;
    p[ccol]= ival;
    if(ic>=iminh){
      // copy all col values
      int pd = ((int*)F->p)[ic];
      int cnz = ((int*)F->nz)[ic];
      nbz[ccol]= cnz;
      for (int i=pd;i<pd+cnz;++i){
        il[ival] = ((int*)F->i)[i];
        xl[ival] = ((double*)F->x)[i];
        ival++;
      }
      ccol++;
    }else{
      break;
    }
    
  }
  
  subfactor * subfact = (subfactor*)malloc( sizeof(subfactor) );
  subfact->size = ccol;
  subfact->i=il;
  subfact->x=xl;
  subfact->nz=nbz;
  subfact->cols=cols;
  subfact->p=p;
  return(subfact);
  
}

void cholmod_tools_restore_subfactor(cholmod_factor * F,subfactor * subF){
  for(int i=0;i<subF->size;++i){
    int istart = subF->p[i];
    int nbz = subF->nz[i];
    int ccol = subF->cols[i];
    //Rcout << "istart :" << istart << ",nbz " << nbz << ", ccol " << ccol << std::endl; 
    ((int *)F->nz)[ccol]=nbz;
    //Rcout << "change nz" << std::endl;
    int istartF = ((int *)F->p)[ccol];
    //Rcout << "istartF : " << istartF << std::endl;
    
    // WARNING CHECK NEEDED
    // ((int *)F->nz)[ccol]>nbz else resize needed ?
    
    for(int iv=0;iv<nbz;++iv){
      ((int*)F->i)[istartF+iv]=subF->i[istart+iv];
      ((double*)F->x)[istartF+iv]=subF->x[istart+iv];
    }
  }
}


double cholmod_tools_cutcost(cholmod_factor * F,double Lintra_diag[],int pivot,int old_pivot,std::set<int> old_pivot_edges,const std::vector<std::pair<int, int>> * cutset,int * permutation,const std::set<int,std::greater<int>> & node_g_intra_nodes,const std::set<int,std::greater<int>> & node_h_intra_nodes, cholmod_common * com){
  
  //double lognbtreeS =  cholmod_tools_logdet_subset(F,node_g_intra_nodes)+cholmod_tools_logdet_subset(F,node_h_intra_nodes);
  std::vector<std::pair<int, int>> cutset_loc(cutset->size());
  int i=0;
  int iminh = old_pivot;
  int iming =  std::numeric_limits<int>::max();
  int ig,ih;
  
  for (auto it=cutset->begin(); it !=cutset->end(); ++it){
    if(node_g_intra_nodes.find(permutation[it->first])!=node_g_intra_nodes.end()){
      ig = permutation[it->first];
      ih = permutation[it->second];
    }else{
      ih = permutation[it->first];
      ig = permutation[it->second];
    }
    if(ig<iming && ig!=pivot){
      iming=ig;
    }
    if(ih<iminh){
      iminh=ih;
    }
    cutset_loc[i] = std::make_pair(permutation[it->first],permutation[it->second]);
    i++;
  }
  
  
  subfactor * subfact = cholmod_tools_save_subfactor(F,node_g_intra_nodes,iming,node_h_intra_nodes,iminh,com);
  //Rcout << "SAVE" << std::endl;
  //print_factor(F);
  
  //double lognbtree=0;
  cholmod_sparse* links1 = cholmod_tools_linkstocol(F->n,pivot,old_pivot,&cutset_loc,com);
  cholmod_updown(true,links1,F,com);


  cholmod_sparse * col1 = cholmod_tools_oldpivotcol(F->n,pivot,old_pivot,old_pivot_edges,&cutset_loc,com);
  cholmod_rowadd(old_pivot,col1,F,com);
  //Rcout << "UP" << std::endl;
  // print_factor(F);
  // Rcout << F->is_monotonic << std::endl;
  // // computelogdet
  double delta_lnbtg = cholmod_tools_logdet_subset_delta_nochange(F,node_g_intra_nodes,Lintra_diag,iming);
  double delta_lnbth = cholmod_tools_logdet_subset_delta_nochange(F,node_h_intra_nodes,Lintra_diag,iminh);
  double delta=delta_lnbtg+delta_lnbth;
  // double lognbtreeT =  cholmod_tools_logdet_subset(F,node_g_intra_nodes)+cholmod_tools_logdet_subset(F,node_h_intra_nodes);
  // Rcout << "----- delta : " << delta << std::endl;
  // Rcout << "----- delta 2 : " << lognbtreeT-lognbtreeS << std::endl;
  // delta= lognbtreeT-lognbtreeS;
  // revert the changes
  cholmod_tools_restore_subfactor(F,subfact);
  //Rcout << "RESTORE" << std::endl;
  //print_factor(F);
  //cholmod_rowdel(old_pivot,NULL,F,com);
  cholmod_free_sparse (&col1, com);
  //cholmod_updown(false,links1,F,com);
  cholmod_free_sparse (&links1, com);
  
  return (delta) ;
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






cholmod_sparse* cholmod_tools_edges_inter(int n, int h, int pivot,std::map<int, int> inter_h,cholmod_common* com ) {
  // links with the pivot ?
  auto search = inter_h.find(pivot);
  double piv_val = 0;
  int link_piv = 0;
  if(search!=inter_h.end()){
    link_piv = 1;
    piv_val = search->second;
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
    if(itr->first>h){
      u=h;
      v=itr->first;
    }else{
      u=itr->first;
      v=h;
    }
    ((int*)S->i)[r*2] = u;
    ((double*)S->x)[r*2] = -sqrt(itr->second);
    ((int*)S->i)[r*2+1] = v;
    ((double*)S->x)[r*2+1] = sqrt(itr->second);
    r++;
  }
  if(link_piv==1){
    ((int*)S->p)[r+1] =((int*)S->p)[r]+1;
    ((int*)S->i)[r*2] = h;
    ((double*)S->x)[r*2] = sqrt(piv_val);
  }
  
  //print_sparse(S);
  return(S);
}


cholmod_sparse* cholmod_tools_edges_diag_inter(int n,const std::map<int, int> & inter_h,cholmod_common* com ) {
  int nb_links = inter_h.size();
  cholmod_sparse *S;
  S = cholmod_allocate_sparse(n,nb_links,nb_links,true,true,0,CHOLMOD_REAL, com);
  int r = 0;
  int u,v;
  ((int*)S->p)[0] = 0;
  for (auto itr = inter_h.begin(); itr!=inter_h.end();++itr){
    ((int*)S->p)[r+1] =((int*)S->p)[r]+1;
    ((int*)S->i)[r] = itr->first;
    ((double*)S->x)[r] = sqrt(itr->second);
    r++;
  }
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

