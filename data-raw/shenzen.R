library(R.matlab)
library(dplyr)
library(sf)
library(ggplot2)
library(gtclust)
nodes_raw  = readMat("./data-raw/NODES.mat")
nodes = data.frame(x=nodes_raw$Node[,1],y=nodes_raw$Node[,2],id=1:nrow(nodes_raw$Node))



links_raw  = readMat("./data-raw/LINKS.mat")
links = data.frame(from=links_raw$Link.1[,1],to=links_raw$Link.1[,2])
links=links |> left_join(nodes,by=c("from"="id"))|> left_join(nodes,by=c("to"="id"),suffix=c("","_to"))
Lp = as.matrix(links[,3:6])
links_geom = st_sfc(lapply(1:nrow(links),\(l){ st_linestring(rbind(Lp[l,1:2],Lp[l,3:4]))}))
links$geometry = links_geom
links.sf=st_sf(links)


speeds_raw  = readMat("./data-raw/SPEEDS.mat")
col_n = paste0("T",floor(((1:288)*5)/60),"H",((1:288)*5)%%60)
days = paste0("D",7:9)
speeds = do.call(cbind,lapply(1:3,\(id){
  df=data.frame(speeds_raw[[id]])
  df[df==0]=NA
  colnames(df)=paste0(days[id],col_n)
  df
}))


speeds_long = speeds |> mutate(link_id=1:n())|> 
  tidyr::pivot_longer(-link_id,names_to = "day_h",values_to = "speed") 
  
  
v=speeds_long |> filter(day_h=="D7T8H30"|day_h=="D7T8H35"|day_h=="D7T8H40") |> 
  group_by(link_id) |> 
  summarise(speed=mean(speed,na.rm=TRUE))  

v$speed[is.na(v$speed)]=mean(v$speed,na.rm = TRUE)
links.sf$speed=v$speed

sum(is.na(links.sf$speed))


ggplot(links.sf)+geom_sf(aes(color=speed*3.6),size=1.1)+
  scale_color_distiller("Speed (km/h)",palette="RdYlGn",direction = 1,limits=c(0,65))+
  theme_void()+
  theme(legend.position="left")+
  ggtitle("Shenzen speed distribution","8h30-8h45")


hc_res=gtclust_lines(links.sf |> select(speed)|>mutate(speed=log(speed)),gtmethod_bayes_dgmm())

gt_res= geocutree(hc_res,9) |> mutate(L=st_length(geometry)) |> arrange(desc(L)) |> head(9)
  
  
ggplot(gt_res)+
  geom_sf(aes(color=factor(cl)))+
  theme_void()+
  scale_color_brewer(palette="Set1",guide="none")+
  ggtitle("Shenzen clustering results","8h30-8h45")

ggplot(gt_res)+
  geom_sf(aes(color=speed*3.6))+
  theme_void()+
  scale_color_distiller("Speed (km/h)",palette="RdYlGn",direction = 1)+
  theme(legend.position="left")+
  ggtitle("Shenzen clustering results","8h30-8h45")


plot(hc_res$Ll)
plot(geocutree(hc_res,4))
plot(hc_res)
nb = sf::st_relate(links.sf,links.sf, pattern = "F***T****")
class(nb)="list"

k_max= 100
N=nrow(links)
pr=sptree_prior(hc_res,k_max)
pr$Cnk=lgamma(N)-lgamma(1:k_max)-lgamma(2013-2:(k_max+1))
pr$ptot=pr$inter+pr$intra-pr$intra[1]#-pr$Cnk-lgamma(1:k_max+1)
pr$Ll=hc_res$Ll[N:(N-k_max+1)]
pr$intra_comp=hc_res$PriorIntra[N:(N-k_max+1)]
pr$inter_comp=hc_res$PriorInter[N:(N-k_max+1)]
pr


library(gtclust)
library(dplyr)
nb=list(c(2,3),c(1,3),c(1,2,4),c(3,5),c(4,6,7),c(5,7),c(5,6))
X=matrix(runif(7*2),nrow=7)

hc_res_small = gtclust_graph(nb,data.frame(X),method = gtmethod_bayes_dgmm(),scaling = "raw")




k_max= 7
N=nrow(X)
pr=sptree_prior(hc_res_small,k_max)
pr$Cnk=lgamma(N)-lgamma(1:k_max)-lgamma(N-2:(k_max+1))
pr$ptot=pr$inter+pr$intra-pr$intra[1]#-pr$Cnk-lgamma(1:k_max+1)
pr$Ll=hc_res_small$Ll[N:(N-k_max+1)]
pr$intra_comp=hc_res_small$PriorIntra[N:(N-k_max+1)]
pr$inter_comp=hc_res_small$PriorInter[N:(N-k_max+1)]
pr |> select(inter,inter_comp,intra,intra_comp)

