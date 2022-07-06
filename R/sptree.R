#' @title Conversion from neighboring list to sparse adjacency matrix
#'
#' @description This function convert a neighboring list to sparse adjacency matrix
#' @param nb neighboring list
#' @return A sparse matrix
#' @export
to_adjmat=function(nb){
  linked = sapply(nb,function(nbs){length(nbs)>0})
  adj_list = do.call(rbind,lapply(which(linked),function(i){data.frame(i=i,j=nb[[i]])}))
  im=max(adj_list)
  Matrix::sparseMatrix(i=adj_list$i,j=adj_list$j,dims = c(im,im))
}

#' @title Compute the laplacian of a graph
#'
#' @description This function compute the laplacian of a graph
#' @param A a sparse adjacency matrix
#' @return A sparse matrix with the laplacian of A
#' @export
laplacian=function(A){
  Matrix::diag(Matrix::rowSums(A))-A
}

#' @title Compute the number of spanning tree in graph A of a graph
#'
#' @description This function compute the number of spanning tree in graph A
#' @param A a sparse adjacency matrix
#' @return the log of the number of spanning tree in graph A
#' @export
log_nb_sptree = function(A){
  if(!is(A,"Matrix")){
    lnbtree = 0;
  }else{
    if(nrow(A)==2){
      if(A[1,2]==0){
        lnbtree = 0
      }else{
        lnbtree = log(A[1,2])
      }
    }else{
      L = laplacian(A)
      lnbtree =   Matrix::determinant(L[-1,-1],log=TRUE)$modulus
    }
  }
  lnbtree
}

decompose_compatible_sptree=function(A,cl){
  pcl = sort(unique(cl))
  K=max(pcl)
  intra_clust_graphs = lapply(pcl,function(k){A[cl==k,cl==k]})
  inter_clust_graph = Matrix::sparseMatrix(i=c(),j=c(),x=1,dims = c(K,K),)
  for (g in 2:K){
    for (h in 1:(g-1)){
        inter_clust_graph[g,h]=sum(A[cl==g,cl==h])
    }
  }
  inter_clust_graph = inter_clust_graph + Matrix::t(inter_clust_graph)
  list(intra_clust_graphs=intra_clust_graphs,inter_clust_graph=inter_clust_graph)
}


log_nb_compatible_sptree = function(A,cl){
  if(all(cl==1)){
    res=list(intra=log_nb_sptree(A),inter=0)
  }else{
    decomposition = decompose_compatible_sptree(A,cl)
    intra = sum(sapply(decomposition$intra_clust_graphs, log_nb_sptree))
    inter = log_nb_sptree(decomposition$inter_clust_graph)
    list(intra = intra,inter=inter)
  }
  
}

#' @title Compute the spanning tree prior of a solution over a range of k value
#'
#' @description This function compute the spanning tree prior of a solution
#' @param sol a gtclust bayesian clustering solution
#' @param k_max the maximum value of k for which the prior should be calculated
#' @return the log of the prior
#' @export
sptree_prior=function(sol,k_max){
  A = to_adjmat(sol$adjacencies_list)
  clp = rep(1,nrow(A))
  res = data.frame(k=1:k_max,inter=NA,intra=NA)
  res[1,"intra"]=log_nb_sptree(A)
  res[1,"inter"]=0
  intra_clust_graphs = list(A)
  inter_clust_graph  = Matrix::sparseMatrix(i=c(),j=c(),x = 0,dims = c(1,1))
  log_nb_intra = c(res[1,"intra"])
  for (g in 2:k_max){
    print(g)
    cl = cutree(sol,g)
    comp_cl = table(clp,cl)
    isplit = which(rowSums(comp_cl>0)==2)
    

    new_intra_clust_graphs_temp=lapply(which(comp_cl[isplit,]>0),function(k){A[cl==k,cl==k]})
    perm = apply(comp_cl,1,function(row){which(row>0)})
    perm_old = unlist(perm[sapply(perm,length)==1])

    
    
    new_clusters=which(comp_cl[isplit,]>0)
    new_intra_clust_graphs=list()
    

    if(g>2){
      new_intra_clust_graphs[perm_old]=intra_clust_graphs[-isplit]
      log_nb_intra[perm_old]=log_nb_intra[-isplit]
    }
    new_intra_clust_graphs[new_clusters]=new_intra_clust_graphs_temp
    intra_clust_graphs = new_intra_clust_graphs
    log_nb_intra[new_clusters] = sapply(intra_clust_graphs[new_clusters], log_nb_sptree)
    res[g,"intra"] = sum(log_nb_intra)
    

    new_inter_clust_graph=Matrix::sparseMatrix(i=c(),j=c(),x = 0,dims = c(g,g))
    if(g>2){
      new_inter_clust_graph[perm_old,perm_old]=inter_clust_graph[-isplit,-isplit]  
    }
    
    for (h in 1:g){
      if(new_clusters[1]!=h){
        ec1 = sum(A[cl==new_clusters[1],cl==h])
        new_inter_clust_graph[new_clusters[1],h] = ec1
        new_inter_clust_graph[h,new_clusters[1]] = ec1
      }
      if(new_clusters[2]!=h){
        ec2 = sum(A[cl==new_clusters[2],cl==h])
        new_inter_clust_graph[new_clusters[2],h] = ec2
        new_inter_clust_graph[h,new_clusters[2]] = ec2
      }

    }
    inter_clust_graph = new_inter_clust_graph
    

    res[g,"inter"]=log_nb_sptree(inter_clust_graph)
    

    clp = cl
  }
  pr=res
  N = length(sol$Ll)
  pr$Cnk=lgamma(N)-lgamma(1:k_max)-lgamma(2013-2:(k_max+1))
  pr$ptot=pr$inter+pr$intra-pr$intra[1]-pr$Cnk-lgamma(1:k_max+1)
  pr$Ll=-sol$Ll[N:(N-k_max+1)]
  pr
}
