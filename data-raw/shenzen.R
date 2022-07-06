library(R.matlab)
library(dplyr)
library(sf)
nodes_raw  = readMat("./data-raw/NODES.mat")
nodes = data.frame(x=nodes_raw$Node[,1],y=nodes_raw$Node[,2],id=1:nrow(nodes_raw$Node))
nodes.sf = 


links_raw  = readMat("./data-raw/LINKS.mat")
links = data.frame(from=links_raw$Link.1[,1],to=links_raw$Link.1[,2])
links=links |> left_join(nodes,by=c("from"="id"))|> left_join(nodes,by=c("to"="id"),suffix=c("","_to"))
Lp = matrix(links[,3:6])
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


ggplot(links.sf)+geom_sf(aes(color=speed))+scale_color_distiller(palette="RdYlGn",direction = 1)+theme_minimal()


hc_res=gtclust_lines(links.sf |> select(speed),gtmethod_bayes_dgmm())
plot(hc_res$Ll)
plot(geocutree(hc_res,7))




