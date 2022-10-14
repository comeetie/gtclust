#ifndef CHOLMOD_TOOLS
#define CHOLMOD_TOOLS

#include "SuiteSparse_config/SuiteSparse_config.h"
#include "CHOLMOD/Include/cholmod.h"
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


void print_sparse(cholmod_sparse * S){
  Rcout << "is packed :" << S->packed << std::endl;
  Rcout << "IX: " ;
  for(int v=0;v<S->nzmax;v++){
    Rcout << "(" << ((int*)S->i)[v] << ", " << ((double*)S->x)[v] << ")";
  }
  Rcout << std::endl;
  Rcout << "P: ";
  for(int v=0;v<(S->ncol+1);v++){
    Rcout << ((int*)S->p)[v] << ", ";
  }
  Rcout << std::endl;
  if(!S->packed){
    Rcout << "NZ: ";
    for(int v=0;v<(S->ncol);v++){
      Rcout << ((int*)S->nz)[v] << ", ";
    }
    Rcout << std::endl;
  }
}


double cholmod_tools_logdet(cholmod_factor * L, cholmod_common * com){
  if(L->is_ll | L->is_super){
    cholmod_change_factor(CHOLMOD_REAL,false,false,true,true,L,com);
  }
  double logdet = 0;
  for(int r=0;r<L->n;r++){

    int pd = ((int*)L->p)[r];
    logdet+= log(((double*)L->x)[pd]);
  }
  return logdet;
}


cholmod_sparse* cholmod_tools_link(int n,int u, int v,   cholmod_common * com) {
  cholmod_sparse *S;
  S = cholmod_allocate_sparse(n,1,2,true,true,0,CHOLMOD_REAL, com);
  ((int*)S->p)[0] = 0;
  ((int*)S->p)[1] = 2;
  ((int*)S->i)[0] = u;
  ((int*)S->i)[1] = v;
  ((double*)S->x)[0] = 1;
  ((double*)S->x)[1] = -1;
  return(S);
}

cholmod_sparse* cholmod_tools_links(int n1,int n2,NumericMatrix cutset,   cholmod_common * com) {
  cholmod_sparse *S;
  
  int nb_links = cutset.nrow();
  S = cholmod_allocate_sparse(n2,nb_links,2*nb_links,true,false,0,CHOLMOD_REAL, com);
  
  ((int*)S->p)[0] = 0;
  int u,v;
  for (int r =0; r< cutset.nrow();r++){
    if(cutset(r,0)>cutset(r,1)){
      u=cutset(r,1);
      v=cutset(r,0);
    }else{
      u=cutset(r,0);
      v=cutset(r,1);
    }
    ((int*)S->p)[r+1] =((int*)S->p)[r]+2;


    // if the link deal with the last node only update u diagonal element
    // ! remember to update n2 with this link
    if(v!=(n2-1) & u!=(n1-1)){
      ((int*)S->i)[r*2+1] = v;
      ((double*)S->x)[r*2+1] = -1;

      ((int*)S->i)[r*2] = u;
      ((double*)S->x)[r*2] = 1;
      
      ((int*)S->nz)[r]=2;
    }else{
      if(v==(n2-1) & u!=(n1-1)){
        ((int*)S->i)[r*2] = u;
        ((double*)S->x)[r*2] = 1;
        ((int*)S->nz)[r]=1;
      }
      
      if(v!=(n2-1) & u==(n1-1)){
        ((int*)S->i)[r*2] = v;
        ((double*)S->x)[r*2] = 1;
        ((int*)S->nz)[r]=1;
      }
      
      if(v==(n2-1) & u==(n1-1)){
        ((int*)S->nz)[r]=0;
      }
      S->packed=false;
    }
  }

  return(S);
}


cholmod_sparse * cholmod_tools_update_col1(int n1,int newn,cholmod_sparse* col, NumericMatrix cutset, cholmod_common * com){
  // deal with specific case newn=2
  if(newn==2){
    cholmod_sparse * new_col = cholmod_allocate_sparse(newn,1,2,true,true,0,CHOLMOD_REAL, com);
    return(new_col);
  }
  int nb_links_incutset = 0;
  int link_with_last = 0;
  int new_links[cutset.nrow()];
  int nl;

  // check the cutset for possible links with the node
  // ! remove links with last node
  for(int r =0; r< cutset.nrow();r++){
    nl = -1;
    if(cutset(r,0)==(n1-1)){
      if(cutset(r,1)!=(newn-1)){
        nl=cutset(r,1);
      }else{
        link_with_last++;
      }
      
    }
    if(cutset(r,1)==(n1-1)){
      if(cutset(r,0)!=(newn-1)){
        nl=cutset(r,0);
      }else{
        link_with_last++;
      }
    }
    if(nl!=-1){
      if(nb_links_incutset==0){
        new_links[0]=nl;
      }else{
        // sorted insertion
        int ii=nb_links_incutset;
        while((new_links[ii]>nl) & (ii>=0)){
          new_links[ii+1]=new_links[ii];
          ii--;
        }
        new_links[ii+1]=nl;
      }
      nb_links_incutset++;
    }
    
  }

  // create new sparse column
  cholmod_sparse * new_col;
  
  if(nb_links_incutset==0 & link_with_last==0){
    // no links in cutset simply copy the old column and shift the rows
    new_col = cholmod_copy_sparse(col,com);
    new_col->nrow=newn;
  }else{
    

    // some links in the cutset merge the two set of links
    new_col = cholmod_allocate_sparse(newn,1,col->nzmax+nb_links_incutset,true,true,0,CHOLMOD_REAL, com);
    // deal with empty col vectors
    if(col->nzmax==1){
      ((double *)new_col->i)[0]=0;
      ((double *)new_col->x)[0]=nb_links_incutset+link_with_last;
    }else{
      for(int v =0; v<(col->nzmax);v++){
        // deal with the diagonal and copy values
        if((((int *)col->i)[v])==(n1-1)){
          ((double *)new_col->x)[v]=((double *)col->x)[v]+nb_links_incutset+link_with_last;
        }else{
          ((double *)new_col->x)[v]=((double *)col->x)[v];
        }
        // copy the rows
        ((int *)new_col->i)[v]=((int *)col->i)[v];
      }
    }    


    
    // insert the newlinks
    for(int v =0; v<nb_links_incutset;v++){
      ((int *)new_col->i)[v+col->nzmax]=new_links[v];
      ((double *)new_col->x)[v+col->nzmax]=-1;
    }
    ((int *)new_col->p)[1]=new_col->nzmax;
  }
  return(new_col);
  
}


cholmod_sparse * cholmod_tools_update_col2(int newn,cholmod_sparse* col, NumericMatrix cutset, cholmod_common * com){
  
  // secial case newn = 2, just two node and one link 
  if(newn==2){
    cholmod_sparse * new_col = cholmod_allocate_sparse(newn,1,2,true,true,0,CHOLMOD_REAL, com);
    ((int*)new_col->p)[0] = 0;
    ((int*)new_col->p)[1] = 2;
    ((int*)new_col->i)[0] = 0;
    ((int*)new_col->i)[1] = 1;
    ((double*)new_col->x)[0] = -1;
    ((double*)new_col->x)[1] = 1;
    return(new_col);
  }
  
  int nb_links_incutset = 0;
  int new_links[cutset.nrow()];
  int nl;
  
  // check the cutset for possible links with the node
  for(int r =0; r< cutset.nrow();r++){
    nl = -1;
    if(cutset(r,0)==newn-1){
      nl=cutset(r,1);
    }
    if(cutset(r,1)==newn-1){
      nl=cutset(r,0);
    }
    if(nl!=-1){
      if(nb_links_incutset==0){
        new_links[0]=nl;
      }else{
        // sorted insertion
        int ii=nb_links_incutset;
        while((new_links[ii]>nl) & (ii>=0)){
          new_links[ii+1]=new_links[ii];
          ii--;
        }
        new_links[ii+1]=nl;
      }
      nb_links_incutset++;
    }
    
  }
  
  // create new sparse column
  cholmod_sparse * new_col;
  int row_shift = newn-col->nrow;
  
  if(nb_links_incutset==0){
    // no links in cutset simply copy the old column and shift the rows
    new_col = cholmod_copy_sparse(col,com);
    new_col->nrow=newn;
    for(int v =0; v<(new_col->nzmax);v++){
      ((int *)new_col->i)[v]=((int *)new_col->i)[v]+row_shift;
    }
  }else{
    // some links in the cutset merge the two set of links
    new_col = cholmod_allocate_sparse(newn,1,col->nzmax+nb_links_incutset,true,true,0,CHOLMOD_REAL, com);
    // insert the newlinks
    for(int v =0; v<nb_links_incutset;v++){
      ((int *)new_col->i)[v]=new_links[v];
      ((double *)new_col->x)[v]=-1;
    }
    
    for(int v =0; v<(col->nzmax);v++){
      // deal with the diagonal and copy values
      if((((int *)col->i)[v]+row_shift)==(newn-1)){
        ((double *)new_col->x)[v+nb_links_incutset]=((double *)col->x)[v]+nb_links_incutset;
      }else{
        ((double *)new_col->x)[v+nb_links_incutset]=((double *)col->x)[v];
      }
      // shift the rows
      ((int *)new_col->i)[v+nb_links_incutset]=((int *)col->i)[v]+row_shift;
    }
    ((int *)new_col->p)[1]=new_col->nzmax;
  }
  return(new_col);
  
}


cholmod_factor * cholmod_tools_combine_factors2(cholmod_factor * L1, cholmod_factor * L2,cholmod_common * com){
  int n = L1->n+L2->n;
  int nzT = L1->nzmax+L2->nzmax;
  cholmod_factor* L = cholmod_allocate_factor(n,com);
  
  L->nzmax=nzT;
  
  void * nz = cholmod_malloc(n, sizeof (int), com) ;
  void * p = cholmod_malloc(n+1, sizeof (int), com) ;
  void * i = cholmod_malloc(nzT, sizeof (int), com) ;
  void * x = cholmod_malloc(nzT, sizeof (double), com) ;
  
  // copy L1  
  size_t sc1 = (L1->n)*sizeof(int);
  memcpy(p,L1->p,sc1);
  memcpy(nz,L1->nz,sc1);
  memcpy(i,L1->i,(L1->nzmax)*sizeof(int));
  memcpy(x,L1->x,(L1->nzmax)*sizeof(double));
  

  // copy L2
  int nbc1 =  L1->n;
  int nbz1 = L1->nzmax;
  for(int c=0;c<=(L2->n);c++){
    ((int*)p)[c+nbc1]=((int*)L2->p)[c]+nbz1;
  }
  memcpy(&((int*)nz)[nbc1],L2->nz,(L2->n)*sizeof(int));
  
  for(int r=0;r<(L2->nzmax);r++){
    ((int*)i)[r+nbz1]=((int*)L2->i)[r]+nbc1;
  }
  memcpy(&((double*)x)[nbz1],L2->x,(L2->nzmax)*sizeof(double));


  // create next/prev array
  void * prev = cholmod_malloc(n+2, sizeof (int), com) ;
  void * next = cholmod_malloc(n+2, sizeof (int), com) ;
  ((int*)prev)[0]=L->n+1;
  ((int*)next)[0]=1;
  for(int c=1;c<L->n;c++){
    ((int*)prev)[c]=c-1;
    ((int*)next)[c]=c+1;
  }
  ((int*)prev)[L->n]=L->n-1;
  ((int*)next)[L->n]=-1;
  
  ((int*)prev)[L->n+1]=-1;
  ((int*)next)[L->n+1]=0;
  

  
  L->p=p;
  L->nz=nz;
  L->x=x;
  L->i=i;
  L->prev=prev;
  L->next=next;
  L->is_ll=false;
  L->is_super=false;
  L->is_monotonic=true;
  L->xtype=CHOLMOD_REAL;
  return L;
}




cholmod_factor * cholmod_tools_Lup(cholmod_factor * F1,cholmod_factor * F2,cholmod_sparse* nr1, std::vector<std::pair<int, int>> cutset,NumericVector permutation, cholmod_common * com){
  

  cholmod_factor *F = cholmod_tools_combine_factors2(F1,F2,com);

  if(F->n>2){
    NumericMatrix cutset_loc(cutset.size(),2);
    int i=0;
    for (auto it=cutset.begin(); it !=cutset.end(); ++it){
      cutset_loc(i,0) = permutation[it->first];
      cutset_loc(i,1) = permutation[it->second];
      i++;
    }

    cholmod_sparse* links1 = cholmod_tools_links(F1->n,F->n,cutset_loc,com);

    cholmod_updown(1,links1,F,com);

    cholmod_free_sparse (&links1, com);
    
    cholmod_sparse * col1 = cholmod_tools_update_col1(F1->n,F->n,nr1,cutset_loc,com);
    cholmod_rowadd(F1->n-1,col1,F,com);

    cholmod_change_factor (CHOLMOD_REAL, F->is_ll, FALSE, TRUE, TRUE, F, com);

    cholmod_free_sparse (&col1, com);
  }

  return (F) ;
}





double cholmod_tools_Lup_logdet(cholmod_factor * F1,cholmod_factor * F2,cholmod_sparse* nr1,NumericMatrix cutset_loc, cholmod_common * com){
  
  
  cholmod_factor *F = cholmod_tools_combine_factors2(F1,F2,com);
  
  if(F->n>2){

    cholmod_sparse* links1 = cholmod_tools_links(F1->n,F->n,cutset_loc,com);
    cholmod_updown(1,links1,F,com);
    cholmod_free_sparse (&links1, com);
    cholmod_sparse * col1 = cholmod_tools_update_col1(F1->n,F->n,nr1,cutset_loc,com);
    cholmod_rowadd(F1->n-1,col1,F,com);
    //cholmod_change_factor (CHOLMOD_REAL, F->is_ll, FALSE, TRUE, TRUE, F, com);
    cholmod_free_sparse (&col1, com);
    
  }
  
  double logdet = cholmod_tools_logdet(F,com);
  cholmod_free_factor(&F, com);
  return (logdet) ;
}



cholmod_sparse * cholmod_tools_colup(int newn,cholmod_sparse* col, std::vector<std::pair<int, int>> cutset,NumericVector permutation, cholmod_common * com){
  NumericMatrix cutset_loc(cutset.size(),2);
  int i=0;
  for (auto it=cutset.begin(); it !=cutset.end(); ++it){
    cutset_loc(i,0) = permutation[it->first];
    cutset_loc(i,1) = permutation[it->second];
    i++;
  }
  return (cholmod_tools_update_col2(newn,col,cutset_loc,com));
  
}

cholmod_factor * init_factor(cholmod_common * com){
  cholmod_factor * L = cholmod_allocate_factor(1,com);
  L->nzmax=1;
  void * nz = cholmod_malloc(1, sizeof (int), com) ;
  void * p = cholmod_malloc(2, sizeof (int), com) ;
  void * i = cholmod_malloc(1, sizeof (int), com) ;
  void * x = cholmod_malloc(1, sizeof (double), com) ;
  
  ((int*)i)[0]=0;
  ((double*)x)[0]=1;
  ((int*)nz)[0]=1;
  ((int*)p)[0]=0;
  ((int*)p)[1]=1;

  L->p=p;
  L->nz=nz;
  L->x=x;
  L->i=i;
  
  void * prev = cholmod_malloc(3, sizeof (int), com) ;
  void * next = cholmod_malloc(3, sizeof (int), com) ;
  ((int*)prev)[0]=L->n+1;
  ((int*)next)[0]=1;

  ((int*)prev)[L->n]=L->n-1;
  ((int*)next)[L->n]=-1;
  
  ((int*)prev)[L->n+1]=-1;
  ((int*)next)[L->n+1]=0;
  
  L->prev=prev;
  L->next=next;

  return(L);
}

#endif


