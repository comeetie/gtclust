library(gtclust)
library(sf)
library(dplyr)
library(ggplot2)


######### FIG 1 (a)

data("modesshare.pts")
df = modesshare.pts |> filter(DEP=="37")

xy = sf::st_coordinates(df)[,1:2]
delaunay = RTriangle::triangulate(RTriangle::pslg(xy))
nb=rep(list(c()),nrow(df))
for (il in 1:nrow(delaunay$E)){
  r = delaunay$E[il,]
  nb[[r[1]]]=c(nb[[r[1]]],r[2])
  nb[[r[2]]]=c(nb[[r[2]]],r[1])
}

Lp = do.call(rbind,lapply(1:length(nb),\(l){
  if(length(nb[[l]])>1){
    cbind(xy[l,1],xy[l,2],xy[nb[[l]],])
  }else{
    c(xy[l,1],xy[l,2],xy[nb[[l]],])
  }
}))

links_geom = st_sfc(lapply(1:nrow(Lp),\(l){ st_linestring(rbind(Lp[l,1:2],Lp[l,3:4]))}),crs=st_crs(modesshare.pts))
nodes_geom = df |> st_centroid() 

xmin = 490000
ymin = 6680000
ggplot()+
  geom_sf(data=st_voronoi(st_union(nodes_geom)),size=0.2,color="#555555",fill="#ffffff")+
  geom_sf(data=links_geom,size=0.2)+
  geom_sf(data=nodes_geom,size=1.5)+
  coord_sf(xlim=c(xmin,xmin+30000),ylim=c(ymin,ymin+21000))+
  theme_void()


######### FIG 1 (b)

data("shenzen")

nb = sf::st_relate(shenzen,shenzen, pattern = "F***T****")
class(nb)="list"
xy = shenzen |> st_centroid() |> st_coordinates()
Lp = do.call(rbind,lapply(1:length(nb),\(l){
  if(length(nb[[l]])>1){
    cbind(xy[l,1],xy[l,2],xy[nb[[l]],])
  }else{
    c(xy[l,1],xy[l,2],xy[nb[[l]],])
  }
}))

links_geom = st_sfc(lapply(1:nrow(Lp),\(l){ st_linestring(rbind(Lp[l,1:2],Lp[l,3:4]))}),crs=st_crs(shenzen))
nodes_geom = shenzen |> st_centroid() 

st_bbox(nodes_geom)

xmin = 2492500
ymin = 38506400
ggplot()+
  geom_sf(data=links_geom,size=0.2)+
  geom_sf(data=nodes_geom,size=2)+
  geom_sf(data=nodes_geom,size=1.5,color="#ffffff")+
  coord_sf(xlim=c(xmin,xmin+2500),ylim=c(ymin,ymin+1750))+
  theme_void()


######### FIG 2 (a)

library(poissoned)
library(Matrix)
library(smoothr)
set.seed(1234)

d=3
v=0.5
ng = 30
nr = 12


centers = rbind(c(d,0),c(-d,0),c(0,d))
xy = do.call(rbind,lapply(1:3,\(g){
  points <- poissoned::poisson_disc(ncols = nr, nrows = nr, cell_size = 10, verbose = TRUE)
  #> poisson_disc(): 500x350, minimum distance = 14.14
  points = t(t(points)-c(nr/2*10,nr/2*10))
  ii=rowSums(points^2)<0.5*max(rowSums(points^2))
  points = points[ii,]/(nr/2*10)*1.5
  points=t(t(points)+centers[g,])
  data.frame(x=points[,1],y=points[,2],g=g)
}))




delaunay = RTriangle::triangulate(RTriangle::pslg(xy[,1:2]))
nb=rep(list(c()),nrow(xy))
for (il in 1:nrow(delaunay$E)){
  r = delaunay$E[il,]
  nb[[r[1]]]=c(nb[[r[1]]],r[2])
  nb[[r[2]]]=c(nb[[r[2]]],r[1])
}

nb

Lp = as.matrix(do.call(rbind,lapply(1:length(nb),\(l){
  if(length(nb[[l]])>1){
    cbind(xy[l,1],xy[l,2],xy[nb[[l]],1:2])
  }else{
    c(xy[l,1],xy[l,2],xy[nb[[l]],1:2])
  }
})))


Lp_ind = do.call(rbind,lapply(1:length(nb),\(l){
  if(length(nb[[l]])>1){
    cbind(l,nb[[l]])
  }else{
    c(l,nb[[l]])
  }
}))



links_geom = st_sfc(lapply(1:nrow(Lp),\(l){ st_linestring(rbind(Lp[l,1:2],Lp[l,3:4]))}))
point_geom = st_sfc(lapply(1:nrow(xy),\(l){st_point(as.numeric(xy[l,1:2]))}))




cl = data.frame(ii=1:nrow(xy),cl=xy[,3])
A=to_adjmat(nb)
ijok = Matrix::which(Matrix::tril(A),arr.ind = TRUE)
cutset = as.data.frame(ijok) |> 
  left_join(cl,by=c("row"="ii")) |> 
  left_join(cl,by=c("col"="ii")) |>
  filter(cl.x!=cl.y)
cc=as.matrix(cutset[,c(1,2)])
A[cc]=FALSE
A[cc[,c(2,1)]]=FALSE
ccsel = cc[c(2,4,9,11),]
A[ccsel]=TRUE
A[ccsel[,c(2,1)]]=TRUE


G=igraph::graph_from_adjacency_matrix(A,mode = "undirected")
ijok = Matrix::which(Matrix::tril(A)==1,arr.ind = TRUE)
ij_undir = data.frame(i=ijok[,1],j=ijok[,2])


links_geom_temp = st_as_sf(data.frame(i=Lp_ind[,1],j=Lp_ind[,2],geom=links_geom)) |> semi_join(ij_undir) 
lnbt_tot = log_nb_sptree(A)
sptree_centrality = rep(0,nrow(Lp))
for (il in 1:nrow(Lp_ind)){
  Ar=A
  Ar[Lp_ind[il,1],Lp_ind[il,2]]=FALSE
  Ar[Lp_ind[il,2],Lp_ind[il,1]]=FALSE
  sptree_centrality[il]=log_nb_sptree(Ar)-lnbt_tot
}


links_geom = st_as_sf(data.frame(i=Lp_ind[,1],j=Lp_ind[,2],centrality=1-exp(sptree_centrality),geom=links_geom)) |> semi_join(ij_undir) 
rc= range(links_geom$centrality)
ggplot()+
  geom_sf(data=links_geom,aes(size=centrality,alpha=centrality))+
  geom_sf(data=point_geom,size=2)+
  geom_sf(data=point_geom,size=1.5,color="#ffffff")+
  scale_alpha_continuous("Spanning tree\ncentrality:",limits=rc)+
  scale_size_continuous("Spanning tree\ncentrality:",limits=rc,range = c(0.5,1.5))+
  theme_void()

ggplot(links_geom|>st_drop_geometry())+
  geom_histogram(aes(x=centrality))+
  labs(title="Spanning tree centrality distribution")+
  theme_bw()


partition_prior=list()
cls=list()
K=3
for(p in 1:1000){
  print(p)
  es=igraph::sample_spanning_tree(G)
  edges_to_rem=sample(length(es),K-1,replace = FALSE)
  forest=igraph::subgraph.edges(G, es,delete.vertices = FALSE)
  forest=igraph::delete.edges(forest,igraph::E(forest)[edges_to_rem])
  co=igraph::components(forest)
  lnbsp=gtclust:::log_nb_compatible_sptree(A,co$membership)
  partition_prior[[p]]=data.frame(intra=lnbsp$intra,inter=lnbsp$inter)  
  cls[[p]]=co$membership
}
prior.df=do.call(rbind,partition_prior) |> 
  mutate(lp=(intra+inter)-lnbt_tot) 



  
pg= st_as_sf(point_geom)
pg$col=factor(cls[[order(prior.df$lp,decreasing = TRUE)[4]]])
pg$col=factor(xy[,3])
ggplot()+
  geom_sf(data=links_geom,aes(size=centrality,alpha=centrality))+
  geom_sf(data=pg,aes(color=col),size=2)+
  scale_alpha_continuous("Spanning tree\ncentrality:",limits=rc)+
  scale_size_continuous("Spanning tree\ncentrality:",limits=rc,range = c(0.5,1.5))+
  theme_void()

ggplot(prior.df)+
  geom_histogram(aes(x=lp))+
  labs(title="Spanning tree centrality distribution")+
  theme_bw()
