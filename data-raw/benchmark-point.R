library(microbenchmark)
library(gtclust)
library(sf)
library(dplyr)
knngraph = function(df,k){
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
  nb
}

delaunaygraph = function(df){
  xy = sf::st_coordinates(df)[,1:2]
  delaunay = RTriangle::triangulate(RTriangle::pslg(xy))
  nb=rep(list(c()),nrow(df))
  for (il in 1:nrow(delaunay$E)){
    r = delaunay$E[il,]
    nb[[r[1]]]=c(nb[[r[1]]],r[2])
    nb[[r[2]]]=c(nb[[r[2]]],r[1])
  }
  nb
}

pts = modesshare.pts |> filter(DEP==75) |> st_union() |> st_geometry() |> st_centroid()
distth = seq(50000,700000,by=50000)
res=list()
for (dt in distth){
  modesshare.pts.sel = st_intersection(modesshare.pts,st_buffer(pts,dt))
  
  time=microbenchmark({
    nb=knngraph(modesshare.pts.sel,5)
  },times = 5)
  res=c(res,list(data.frame(dt=dt,
                            N=nrow(modesshare.pts.sel),
                            M=sum(sapply(nb, length)),
                            alg="knn",
                            time=median(time$time))))
  
  time=microbenchmark({
    hc=gtclust_graph(nb,modesshare.pts.sel)
  },times = 5)
  res=c(res,list(data.frame(dt=dt,
                            N=nrow(modesshare.pts.sel),
                            M=sum(sapply(nb, length)),
                            alg="gtclust-knn",
                            time=median(time$time))))
  
  
  time=microbenchmark({
    hc=gtclust_graph(nb,modesshare.pts.sel,method="bayes_mom")
  },times = 5)
  res=c(res,list(data.frame(dt=dt,
                            N=nrow(modesshare.pts.sel),
                            M=sum(sapply(nb, length)),
                            alg="gtclust-knn-bayes",
                            time=median(time$time))))
  
  time=microbenchmark({
    nb=delaunaygraph(modesshare.pts.sel)
  },times = 5)
  res=c(res,list(data.frame(dt=dt,
                            N=nrow(modesshare.pts.sel),
                            M=sum(sapply(nb, length)),
                            alg="delaunay",
                            time=median(time$time))))
  
  time=microbenchmark({
    hc=gtclust_graph(nb,modesshare.pts.sel)
  },times = 5)
  res=c(res,list(data.frame(dt=dt,
                            N=nrow(modesshare.pts.sel),
                            M=sum(sapply(nb, length)),
                            alg="gtclust-delaunay",
                            time=median(time$time))))
  
  time=microbenchmark({
    hc=gtclust_graph(nb,modesshare.pts.sel,method="bayes_mom")
  },times = 5)
  res=c(res,list(data.frame(dt=dt,
                            N=nrow(modesshare.pts.sel),
                            M=sum(sapply(nb, length)),
                            alg="gtclust-delaunay-bayes",
                            time=median(time$time))))
  
}
res.bench = do.call(rbind,res)