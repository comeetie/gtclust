#' gtclust: A package for fast clustering of spatial or temporal data with contiguity constrained hierarchical clustering
#'
#' The mypackage package provides three categories of important functions:
#' gtclust_graph, gtclust_temp, gt_poly.
#' 
#' @section gtclust functions:
#' The mypackage functions ...
#'
#' @docType package
#' @name gtclust
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib gtclust, .registration=TRUE
NULL
#> NULL


#' @title Hierarchical clustering with contiguity constraints for temporal data
#'
#' @description This function take a data.frame and performs hierarchical clustering with contiguity constraints.
#' @param df an \code{\link[sf]{sf}} data.frame with polygons like features
#' @param method linkage criterion in ward (default) or average, median
#' @param scaling default scaling of the features in zscore (default) or raw (i.e. no scaling)
#' @return an \code{\link[stats]{hclust}} like object with additional slots
#' \describe{
#'   \item{data}{The numeric data (eventually scaled) used for the clustering}
#'   \item{centers}{The protoypes of each tree nodes}
#' }
#' @export
gtclust_temp=function(df,method="ward",scaling="raw"){
  nT = nrow(df)
  nb = lapply(1:nT,\(it){
    nei = c(it-1,it+1)
    nei[nei>0 & nei<=nT]
  })
  hc_res=gtclust_graph(nb,df,method,scaling)
  hc_res$call=sys.call()
  hc_res
}



#' @title Hierarchical clustering with contiguity constraints for point data with delaunay links 
#'
#' @description This function take a data.frame and performs hierarchical clustering with contiguity constraints.
#' @param df an \code{\link[sf]{sf}} data.frame with polygons like features
#' @param method linkage criterion in ward (default) or average, median
#' @param scaling default scaling of the features in zscore (default) or raw (i.e. no scaling)
#' @return an \code{\link[stats]{hclust}} like object with additional slots
#' \describe{
#'   \item{data}{The numeric data (eventually scaled) used for the clustering}
#'   \item{centers}{The protoypes of each tree nodes}
#'   \item{leafs_geometry}{geometries of the dendrogram leafs as an sfc list}
#'   \item{geotree}{geometries of the dendrogram no-leafs node as an sfc list}
#' }
#' @export
gtclust_delaunay=function(df,method="ward",scaling="raw"){
  if(!methods::is(df,"sf")){
    stop("The dataset must be an sf data.frame.",call. = FALSE)
  }
  
  if(!all(sapply(sf::st_geometry(df),function(u){sf::st_is(u,"POINT")}))){
    stop("The dataset must contains only POINTS.",call. = FALSE)
  }
  df_nogeo=sf::st_drop_geometry(df)
  
  xy = sf::st_coordinates(df)[,1:2]
  delaunay = RTriangle::triangulate(RTriangle::pslg(xy))
  nb=rep(list(c()),nrow(df))
  for (il in 1:nrow(delaunay$E)){
    r = delaunay$E[il,]
    nb[[r[1]]]=c(nb[[r[1]]],r[2])
    nb[[r[2]]]=c(nb[[r[2]]],r[1])
  }

  hc_res=gtclust_graph(nb,df_nogeo,method,scaling)
  hc_res$call=sys.call()
  # add geographical data
  hc_res$leafs_geometry = sf::st_geometry(df)
  hc_res$geotree = build_geotree(hc_res$merge,df)
  class(hc_res)=c(class(hc_res),"geoclust")
  hc_res
}


#' @title Hierarchical clustering with contiguity constraints for point data with distance threshold links 
#'
#' @description This function take a data.frame and performs hierarchical clustering with contiguity constraints.
#' @param df an \code{\link[sf]{sf}} data.frame with polygons like features
#' @param epsilon maximum distance allowed
#' @param method linkage criterion in ward (default) or average, median
#' @param scaling default scaling of the features in zscore (default) or raw (i.e. no scaling)
#' @return an \code{\link[stats]{hclust}} like object with additional slots
#' \describe{
#'   \item{data}{The numeric data (eventually scaled) used for the clustering}
#'   \item{centers}{The protoypes of each tree nodes}
#'   \item{leafs_geometry}{geometries of the dendrogram leafs as an sfc list}
#'   \item{geotree}{geometries of the dendrogram no-leafs node as an sfc list}
#' }
#' @export
gtclust_dist=function(df,epsilon,method="ward",scaling="raw"){
  if(!methods::is(df,"sf")){
    stop("The dataset must be an sf data.frame.",call. = FALSE)
  }
  
  if(!all(sapply(sf::st_geometry(df),function(u){sf::st_is(u,"POINT")}))){
    stop("The dataset must contains only POINTS.",call. = FALSE)
  }
  df_nogeo=sf::st_drop_geometry(df)
  
  # buid graph
  buf = sf::st_buffer(df,epsilon)
  nb = sf::st_intersects(df,buf)
  class(nb)="list"
  
  hc_res=gtclust_graph(nb,df_nogeo,method,scaling)
  hc_res$call=sys.call()
  # add geographical data
  hc_res$leafs_geometry = sf::st_geometry(df)
  hc_res$geotree = build_geotree(hc_res$merge,df)
  class(hc_res)=c(class(hc_res),"geoclust")
  hc_res
}

#' @title Hierarchical clustering with contiguity constraints for point data with knn links 
#'
#' @description This function take a data.frame and performs hierarchical clustering with contiguity constraints.
#' @param df an \code{\link[sf]{sf}} data.frame with polygons like features
#' @param k number of nearest neighbors to take for building the graph (the graph will be symmetric so some points may have in fine more neighbors)
#' @param method linkage criterion in ward (default) or average, median
#' @param scaling default scaling of the features in zscore (default) or raw (i.e. no scaling)
#' @return an \code{\link[stats]{hclust}} like object with additional slots
#' \describe{
#'   \item{data}{The numeric data (eventually scaled) used for the clustering}
#'   \item{centers}{The protoypes of each tree nodes}
#'   \item{leafs_geometry}{geometries of the dendrogram leafs as an sfc list}
#'   \item{geotree}{geometries of the dendrogram no-leafs node as an sfc list}
#' }
#' @export
gtclust_knn=function(df,k=3,method="ward",scaling="raw"){
  if(!methods::is(df,"sf")){
    stop("The dataset must be an sf data.frame.",call. = FALSE)
  }
  
  if(!all(sapply(sf::st_geometry(df),function(u){sf::st_is(u,"POINT")}))){
    stop("The dataset must contains only POINTS.",call. = FALSE)
  }
  df_nogeo=sf::st_drop_geometry(df)
  
  # build graph
  xy = sf::st_coordinates(df)[,1:2]  
  knn = RANN::nn2(xy,k=k)
  # ensure symmetry and extract adjacency list from results
  nb = rep(list(c()),nrow(df))
  for (i in 1:nrow(xy)){
    knei = setdiff(knn$nn.idx[i,],i)
    nb[[i]]=unique(c(nb[[i]],knei))
    for(j in knn$nn.idx[i,]){
      nb[[j]]=unique(c(nb[[j]],i))
    }
  }

  
  hc_res=gtclust_graph(nb,df_nogeo,method,scaling)
  hc_res$call=sys.call()
  # add geographical data
  hc_res$leafs_geometry = sf::st_geometry(df)
  hc_res$geotree = build_geotree(hc_res$merge,df)
  class(hc_res)=c(class(hc_res),"geoclust")
  hc_res
}



#' @title Hierarchical clustering with contiguity constraints between polygons
#'
#' @description This function take an \code{\link[sf]{sf}} data.frame and performs hierarchical clustering with contiguity constraints.
#' @param df an \code{\link[sf]{sf}} data.frame with polygons like features
#' @param method linkage criterion in ward (default) or average, median
#' @param scaling default scaling of the features in zscore (default) or raw (i.e. no scaling)
#' @param adjacency adjacency type to use  "rook" (default) or queen
#' @return an \code{\link[stats]{hclust}} like object with additional slots
#' \describe{
#'   \item{leafs_geometry}{geometries of the dendrogram leafs as an sfc list}
#'   \item{geotree}{geometries of the dendrogram no-leafs node as an sfc list}
#'   \item{data}{The numeric data (eventually scaled) used for the clustering}
#'   \item{centers}{The protoypes of each tree nodes}
#' }
#' @export
gtclust_poly=function(df,method="ward",adjacency="rook",scaling="raw"){
  
  if(!methods::is(df,"sf")){
    stop("The dataset must be an sf data.frame.",call. = FALSE)
  }
  
  if(!all(sapply(sf::st_geometry(df),function(u){sf::st_is(u,"MULTIPOLYGON")}) | 
          sapply(sf::st_geometry(df),function(u){sf::st_is(u,"POLYGON")}))){
    stop("The dataset must contains only POLYGONS or MULTIPOLYGONS.",call. = FALSE)
  }
  df_nogeo=sf::st_drop_geometry(df)

  
  # build graph
  # see https://github.com/r-spatial/sf/issues/234#issuecomment-300511129 and ?st_relate
  if(adjacency=="rook"){
    nb = sf::st_relate(df, df, pattern = "F***1****")
  }else{
    nb = sf::st_relate(df,df, pattern = "F***T****")
  }
  class(nb)="list"
  hc_res=gtclust_graph(nb,df_nogeo,method,scaling)
  hc_res$call=sys.call()
  # add geographical data
  hc_res$leafs_geometry = sf::st_geometry(df)
  hc_res$geotree = build_geotree(hc_res$merge,df)
  class(hc_res)=c(class(hc_res),"geoclust")
  hc_res
}




#' @title Hierarchical clustering with contiguity constraints
#'
#' @description This function take an data.frame and performs hierarchical clustering with contiguity 
#' constraints using a graph describing the contiguity (provided )
#' @param adjacencies_list graph describing the contiguity between the rows of df as a list of adjacencies 
#' @param df a data.frame with numeric columns
#' @param method linkage criterion in ward (default) or average, median
#' @param scaling default scaling of the features in zscore or raw (i.e. no scaling, the default)
#' @return an \code{\link[stats]{hclust}} like object with additional slots
#' \describe{
#'   \item{data}{The numeric data (eventually scaled) used for the clustering}
#'   \item{centers}{The protoypes of each tree nodes}
#' }
#' @export
gtclust_graph = function(adjacencies_list,df,method="ward",scaling="raw"){
  
  
  if(is.character(method) && !(method %in% c("ward","centroid","median","chisq","bayes_mom","bayes_dgmm"))){
    stop("'method' not recognized")
  }
  
  
  if (is.character(method)) 
    method <- get(paste0("gtmethod_",method), mode = "function", envir = parent.frame())
  if (is.function(method)) 
    method <- method()
  if (is.null(method$method)) {
    print(method)
    stop("'method' not recognized")
  }
  
  if(!(scaling %in% c("zscore","raw"))){
    stop("The scaling argument must be zscore or raw.")
  }
  if(!(methods::is(df,"data.frame")|methods::is(df,"matrix"))){
    stop("df must be a data.frame or a matrix")
  }
  if(methods::is(df,"matrix") & !is.numeric(df)){
    stop("df must be numeric.")
  }
  
  if(methods::is(df,"data.frame")){
    # remove geo in case
    if(methods::is(df,"sf")){
      df= sf::st_drop_geometry(df)
    }
    # select only numeric features
    num_feats = unlist(lapply(df,is.numeric))
    if(sum(num_feats)!=ncol(df)){
      warning("Some features were not numeric and have been removed from the clustering.",call. = FALSE)
      df=df[,num_feats]
    }
  }

  # check for missing values
  if(sum(is.na(df))>0){
    stop("Some regions have missing values and missing values are not allowed.",call. = FALSE)
  }
  
  # scales
  if(scaling=="zscore"){
    df_scaled = apply(df,2,\(col){(col-mean(col))/stats::sd(col)})
  }else{
    df_scaled = as.matrix(df)
  }
  
  nb_c = lapply(adjacencies_list,\(nei){nei-1})
  # run the algorithm
  res=hclustcc_cpp(nb_c,df_scaled,method)

  

  # format the results in hclust form
  hc_res = list(merge=res$merge, 
                height=compute_height(cumsum(res$height)),
                order=order_tree(res$merge,nrow(res$merge)),
                labels=(rownames(df)),
                call=sys.call(),
                method=method$method,
                dist.method="euclidean",
                k.relaxed=res$k.relaxed,
                data=res$data,
                centers=res$centers)
  if(methods::is(method,"bayesian_gtmethod")){
    hc_res = c(hc_res,list(Ll=res$Ll))
  }
  class(hc_res)  <- c("gtclust","hclust")

  hc_res
}


#' @title Cut a Geograpĥic Tree into Groups of Data and return an sf data.frame 
#'
#' @description Cuts a tree, e.g., as resulting from geohclust_poly, into several groups either by specifying the desired number(s) of groups or the cut height(s).
#' @param tree a tree as produced by gtclust. cutree() only expects a list with components merge, height, and labels, of appropriate content each.
#' @param k an integer scalar or vector with the desired number of groups
#' @param h numeric scalar or vector with heights where the tree should be cut.
#' At least one of k or h must be specified, k overrides h if both are given.
#' @return an \code{\link[sf]{sf}} like object
#' @export
geocutree=function(tree,k = NULL, h= NULL){
  if(!methods::is(tree,"geoclust")){
    stop("geocutree only accepts gtclust objects.")
  }
  if(is.null(k) && is.null(h)){
    stop("At least one of k or h must be specified, k overrides h if both are given.")
  }
  N=nrow(tree$merge)+1
  if(!is.null(h) & is.null(k)){
    k <- N + 1L - apply(outer(c(tree$height, Inf), h, `>`),2, which.max)
  }
  cl = stats::cutree(tree,k=k)

  istart = which(!duplicated(cl))
  clust_geo = list()
  clust_x = matrix(0,nrow=k,ncol=ncol(tree$data))
  ck=1
  for (i in istart){
    f = -i
    cnode = -i
    while(length(f)!=0){
      f = which(tree$merge[1:(N-k),1]==f | tree$merge[1:(N-k),2]==f)
      if(length(f)>0){
        cnode = f
      }
    }
    if(cnode>0){
      clust_geo[[ck]]=tree$geotree[[cnode]]
      clust_x[ck,]=tree$centers[cnode,]
    }else{
      clust_geo[[ck]]=tree$leafs_geometry[[-cnode]]
      clust_x[ck,]=tree$data[-cnode,]
    }
    ck=ck+1
  }
  colnames(clust_x)=colnames(tree$centers)
  sf::st_sf(cl=1:k,clust_x,geometry=sf::st_as_sfc(clust_geo,crs=sf::st_crs(tree$leafs_geometry)))
}




order_tree=function(merge,i){
  if(merge[i,1]<0){
    left = -merge[i,1];
  }else{
    left = order_tree(merge,merge[i,1])
  }
  if(merge[i,2]<0){
    right = -merge[i,2];
  }else{
    right = order_tree(merge,merge[i,2])
  }
  c(left,right)
}


build_geotree=function(merge,df){
  geoms= sf::st_geometry(df)
  geotree = list()
  for (i in 1:nrow(merge)){
    if(merge[i,1]<0){
      left = geoms[[-merge[i,1]]];
    }else{
      left = geotree[[merge[i,1]]]
    }
    if(merge[i,2]<0){
      right = geoms[[-merge[i,2]]];
    }else{
      right = geotree[[merge[i,2]]]
    }
    geotree[[i]] = sf::st_union(left,right)
  }
  geotree
}



#' @title aggregation methods for gtclust
#' @return An object of class gtmethod
#' @describeIn gtmethod classical ward method
#' @export
gtmethod_ward = function(){
  structure(list(method = "ward"), 
            class = "gtmethod")
}

#' @describeIn gtmethod classical centroid method
#' @export
gtmethod_centroid = function(){
  structure(list(method = "centroid"), 
            class = "gtmethod")
}


#' @describeIn gtmethod classical median method
#' @export
gtmethod_median = function(){
  structure(list(method = "median"), 
            class = "gtmethod")
}

#' @describeIn gtmethod chi-square method
#' @export
gtmethod_chisq = function(){
  structure(list(method = "chisq"), 
            class = "gtmethod")
}


#' @param beta prior parameters for the dirichlet distribution 
#' @describeIn gtmethod bayesian mixture of multinomials
#' @export
gtmethod_bayes_mom = function(beta = 1){
  structure(list(method = "bayes_mom",beta = beta), 
            class = c("gtmethod","bayesian_gtmethod"))
}

#' @param tau Prior parameter (inverse variance), (default 0.01)
#' @param kappa Prior parameter (gamma shape), (default to 1)
#' @param beta Prior parameter (gamma rate), (default to NaN, in this case beta will be estimated from data as 0.1 time the mean of X columns variances)
#' @param mu Prior for the means (vector of size D), (default to NaN, in this case mu will be estimated from data as the mean of X)
#' @describeIn gtmethod bayesian diagonal gaussian mixture model
#' @export
gtmethod_bayes_dgmm = function(tau = 0.01, kappa = 1, beta = NaN, mu = as.matrix(NaN)){
  structure(list(method = "bayes_dgmm",tau=tau,kappa=kappa,beta = beta,mu=mu), 
            class = c("gtmethod","bayesian_gtmethod"))
}

#' @title plot gtclust dendrogram
#'
#' @param x a tree as produced by a gtclust variant. 
#' @param nb_max_leafs number of leafs to keep 
#' @return an \code{\link[ggplot2]{ggplot}} object
#' @export
plot.gtclust=function(x,y=NULL,nb_max_leafs=500,...){
  
  tree = x 
  rownames(tree$data)=1:nrow(tree$data)
  if(substr(tree$method,1,5)=="bayes"){
    im = which.min(tree$Ll)
    nb_max_leafs=(length(tree$height)+1)-im
  }else{
    im=(length(tree$height)+1)-nb_max_leafs
  }
  small_tree = maptree::clip.clust(tree,data=tree$data,k = nb_max_leafs)

  dend_data = ggdendro::dendro_data(as.dendrogram(small_tree))
  cluster_sizes = sapply(small_tree$membership,length)
  dend_data$labels$size=cluster_sizes
  
  if(tree$k.relaxed>1){
    ymax =tree$height[nrow(tree$data)-tree$k.relaxed]
  }else{
    ymax=max(small_tree$height)
  }

  segs_noconstr = dend_data$segments[dend_data$segments$y>ymax,]
  segs_constr = dend_data$segments[dend_data$segments$y<=ymax,]
  ggplot2::ggplot() + 
    ggplot2::geom_segment(data=segs_constr,ggplot2::aes_(x =~ x, y  =~ y, xend =~ xend, yend =~ yend,linetype="constrained"),size=0.3)+
    ggplot2::geom_segment(data=segs_noconstr,ggplot2::aes_(x =~ x, y =~ y, xend =~ xend, yend =~ yend,linetype="relaxed merge"),color="#aeaeae",size=0.5,linetype="dotted")+
    ggplot2::geom_point(data = dend_data$labels, ggplot2::aes_(x=~x, y=~y-0.01*max(small_tree$height), size=~size))+
    ggplot2::theme_bw()+
    ggplot2::scale_x_continuous("",breaks=c())+
    ggplot2::scale_linetype_manual(values=c("relaxed"="dotted","constrained"="solid"))+
    ggplot2::scale_y_continuous("Height",n.breaks = 8)+
    ggplot2::guides(linetype=ggplot2::guide_legend("Merge type:"),size=ggplot2::guide_legend("Branch size:"))+
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      legend.position="bottom"
    )
  
  
}

collapse = function(tree,nb_max_leaf=500){
  L=length(tree$height)
  cl=cutree(tree,nb_max_leaf)
  itokeep = seq(L-nb_max_leaf+2,L)
  merge=tree$merge[itokeep,]
  nodes=as.vector(merge)
  leafs = !(nodes %in% itokeep)
  sum(leafs)
  nodes[leafs]=-seq(1,sum(leafs))
  nodes[!leafs]=nodes[!leafs]-min(nodes[!leafs])+1
  merge=matrix(nodes,ncol=2)
  tree$merge=merge
  tree$order=gtclust:::order_tree(tree$merge,1)
  tree$height=tree$height[itokeep]
  tree$members=lapply(seq_len(nb_max_leaf),function(k){which(cl==k)})
  tree
}




# extract the pareto front
compute_height <- function(heights) {
  # vector of icls value from root to leaves
  icl <- -heights[length(heights):1]
  # K
  K <- 1:length(icl)
  
  # init H
  H <- rep(0, length(icl))
  
  # current merge position
  cdi <- Inf
  # current best line
  bestline <- 1
  # vector with indexes of solutions that belong to the pareto front
  Front <- c(1)
  
  # from root to leaves
  for (l in 2:length(icl)) {
    
    # merge value with current bestline
    di <- (icl[l] - icl[bestline])
    din <- di / (l - bestline)
    
    # is their a potential merge ?
    if (di > 0) {
      
      # if this merge did not occurs after the current one update the front
      while (din > cdi & length(Front) > 1) {
        
        # remove the last solution from the front
        Front <- Front[-length(Front)]
        H[bestline] <- -1
        
        # update bestline
        bestline <- Front[length(Front)]
        
        # update merge position
        di <- (icl[l] - icl[bestline])
        din <- di / (l - bestline)
        # update previous merge position
        if (length(Front) > 1) {
          cdi <- (icl[bestline] - icl[Front[length(Front) - 1]]) / (bestline - Front[length(Front) - 1])
        } else {
          cdi <- H[1]
        }
      }
      
      # add the extracted solution to the front
      H[Front[length(Front)]] <- din
      cdi <- din
      bestline <- l
      Front <- c(Front, l)
    } else {
      # if solution not in front
      H[l] <- -1
    }
  }
  
  # copy from left previous value
  for (l in 2:length(icl)) {
    if (H[l] == -1) {
      H[l] <- H[l - 1]
    }
  }
  H[length(H):1]
}
