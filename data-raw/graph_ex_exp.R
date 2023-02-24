library(gtclust)
library(sf)
library(sfnetworks)
library(dplyr)
library(stringr)
library(ggplot2)


######### FIG 1 (a)

data("modesshare.rc")
df = modesshare.rc |> filter(str_sub(INSEE_COM,1,2)=="37")

xy = sf::st_coordinates(df|> st_centroid())[,1:2]
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
  geom_sf(data=df,size=0.2,color="#555555",fill="#ffffff")+
  geom_sf(data=links_geom,size=0.2)+
  geom_sf(data=nodes_geom,size=2)+
  geom_sf(data=nodes_geom,size=1.5,color="white")+
  coord_sf(xlim=c(xmin,xmin+30000),ylim=c(ymin,ymin+21000))+
  theme_void()
ggsave("./data-raw/images/graph_ex_1.pdf",width=5,height=4)


######### FIG 1 (b)

data("shenzhen.net")
shenzhen = st_as_sf(shenzhen.net,"edges")
nb = sf::st_relate(shenzhen,shenzhen, pattern = "F***T****")
class(nb)="list"
xy = shenzhen |> st_line_sample(1) |> st_cast("POINT") |> st_coordinates()
Lp = do.call(rbind,lapply(1:length(nb),\(l){
  if(length(nb[[l]])>1){
    cbind(xy[l,1],xy[l,2],xy[nb[[l]],])
  }else{
    c(xy[l,1],xy[l,2],xy[nb[[l]],])
  }
}))

links_geom = st_sfc(lapply(1:nrow(Lp),\(l){ st_linestring(rbind(Lp[l,1:2],Lp[l,3:4]))}),crs=st_crs(shenzhen))
nodes_geom = shenzhen |> st_line_sample(1) |> st_cast("POINT")

st_bbox(nodes_geom)

xmin = 2492500
ymin = 38506400
ggplot()+
  geom_sf(data=st_as_sf(shenzhen.net,"edges"),size=0.4)+
  geom_sf(data=nodes_geom,size=2)+
  geom_sf(data=nodes_geom,size=1.5,color="#ffffff")+
  geom_sf(data=st_as_sf(shenzhen.net,"nodes"),size=2)+
  geom_sf(data=st_as_sf(shenzhen.net,"nodes"),size=1.5,color="#cc2222")+
  coord_sf(xlim=c(xmin,xmin+2500),ylim=c(ymin,ymin+1750))+
  theme_void()
ggsave("./data-raw/images/graph_ex_2.pdf",width=5,height=4)