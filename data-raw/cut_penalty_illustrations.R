library(osmdata)

library(igraph)
library(gtclust)
library(ggplot2)
library(sf)
library(dplyr)
# 
# roads <- opq(bbox = 'Paris fr') %>%
#   add_osm_feature(key = 'highway') %>%
#   osmdata_sf()
# plot(roads$osm_lines |> st_geometry())
# 
# roads_amboise = opq(bbox = 'Amboise fr') %>%
#   add_osm_feature(key = 'highway') %>%
#   osmdata_sf()

plot(roads_amboise$osm_lines |> st_geometry())

roads_amboise_small = read_sf("./data-raw/roads_amboise_small.gpkg")

nb = sf::st_intersects(roads_amboise_small,roads_amboise_small)
class(nb)="list"


A=to_adjmat(nb)
G=igraph::graph_from_adjacency_matrix(A)
co = components(G)

roads_connected = roads_amboise_small[co$membership==1,]


ggplot(roads_connected)+geom_sf(size=0.3)+theme_void()

nb = sf::st_intersects(roads_connected,roads_connected)
class(nb)="list"
A=to_adjmat(nb)
G=igraph::graph_from_adjacency_matrix(A)


roads_amboise_as_node = roads_connected |> st_centroid()

#roads_amboise_as_node = roads_connected |> st_cast("POINT") |> filter(!duplicated(osm_id))

XY=st_coordinates(roads_amboise_as_node)[,1:2]


nodes=roads_amboise_as_node |> select(osm_id,name)

Lp = do.call(rbind,lapply(1:length(nb),\(l){
  cbind(XY[l,1],XY[l,2],XY[nb[[l]],])
}))

links_geom = st_sfc(lapply(1:nrow(Lp),\(l){ st_linestring(rbind(Lp[l,1:2],Lp[l,3:4]))}),crs=4326)



ggplot(links_geom)+geom_sf(size=0.3)+geom_sf(data=roads_amboise_as_node,size=0.3)+theme_void()


spT=minimum.spanning.tree(G)
nbe= ecount(spT)
spTcut=delete_edges(spT,edges = sample(nbe,1))
co=components(spTcut)
while(co$csize[2]<500){
  spTcut=delete_edges(spT,edges = sample(nbe,1))
  co=components(spTcut)
}
co
nodes$cl = co$membership

lnbt_tot = log_nb_sptree(A)
A1=A[co$membership==1,co$membership==1]
lnbt_A1 = log_nb_sptree(A1)
A2=A[co$membership==2,co$membership==2]
lnbt_A2 = log_nb_sptree(A2)
l_ncutset=log(sum(A[co$membership==1,co$membership==2]))

cut_cost = lnbt_tot-lnbt_A1-lnbt_A2-l_ncutset

cut_cost
exp(cut_cost)

ggplot(links_geom|>st_transform(2154))+
  geom_sf(size=0.3,color="#bbbbbb")+
  geom_sf(data=nodes|>st_transform(2154),aes(color=factor(cl)),size=1)+
  theme_void()+
  scale_color_brewer(palette = "Set1",guide="none")+
  ggtitle(paste0("Cut cost: ",round(cut_cost,2)),paste0("so as many as ",round(exp(cut_cost)/1000),"(K) times less spanning trees compatible with the cut than without"))


