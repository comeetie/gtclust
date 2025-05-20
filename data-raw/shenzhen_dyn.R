library(R.matlab)
library(dplyr)
library(sf)
library(ggplot2)
library(sfnetworks)
library(gtclust)
library(igraph)
library(stringr)
library(tidygraph)
library(lubridate)
library(ggpubr)
# see https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0260201
# and also https://figshare.com/articles/dataset/Shenzhen_whole_day_Speeds/7212230

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
  tidyr::pivot_longer(-link_id,names_to = "day_h",values_to = "speed") |>
  mutate(day=str_sub(day_h,1,2),hm=str_sub(day_h,3,20)) |>
  tidyr::separate(hm,sep = "H",into = c("H","M")) |>
  mutate(H=as.numeric(gsub("T","",H)),M=as.numeric(M))
  


v_pick = speeds_long |> filter(H>=5,H<=8) |>   
  group_by(day,link_id) |> 
  summarise(speed=mean(speed,na.rm=TRUE)) |>
  tidyr::pivot_wider(names_from = day,values_from = speed) |> arrange(link_id)

v_pick


shenhzen_pick=bind_cols(links.sf,v_pick)
shenzhen.net.pick = as_sfnetwork(shenhzen_pick)
shenzhen.net.pick = convert(shenzhen.net.pick, to_spatial_smooth,summarise_attributes=list(D7="mean",D8="mean",D9="mean"))
usethis::use_data(shenzhen.net.pick,overwrite = TRUE)



hc=gtclust_lines(st_as_sf(shenzhen.net.pick,"edges")|>select(D7)|>mutate(D7=3.6*D7),method = gtmethod_bayes_dgmm(beta = 15,tau=0.0001))
d7=c()
for (K in 2:5){
  df = data.frame(v=shenzhen.net.pick |> activate("edges") |> pull(D8),cl=cutree(hc,K))
  tv=df |> group_by(cl) |> summarise(vs=var(v),n=n()) |> mutate(tv=vs*n) |> summarise(tv=sum(tv)) /(nrow(df)*var(df$v))
  d7=c(d7,tv$tv)
}


hc=gtclust_lines(st_as_sf(shenzhen.net.pick,"edges")|>select(D8)|>mutate(D8=3.6*D8),method = gtmethod_bayes_dgmm(beta = 15,tau=0.0001))
d8=c()
for (K in 2:5){
  df = data.frame(v=shenzhen.net.pick |> activate("edges") |> pull(D8),cl=cutree(hc,K))
  tv=df |> group_by(cl) |> summarise(vs=var(v),n=n()) |> mutate(tv=vs*n) |> summarise(tv=sum(tv)) /(nrow(df)*var(df$v))
  d8=c(d8,tv$tv)
}

hc=gtclust_lines(st_as_sf(shenzhen.net.pick,"edges")|>select(D9)|>mutate(D9=3.6*D9),method = gtmethod_bayes_dgmm(beta = 15,tau=0.0001))
d9=c()
for (K in 2:5){
  df = data.frame(v=shenzhen.net.pick |> activate("edges") |> pull(D8),cl=cutree(hc,K))
  tv=df |> group_by(cl) |> summarise(vs=var(v),n=n()) |> mutate(tv=vs*n) |> summarise(tv=sum(tv)) /(nrow(df)*var(df$v))
  d9=c(d9,tv$tv)
}

data.frame(K=2:5,d7=d7,d8=d8,d9=d9)



###############################
### Clustering temp 
speed_quarter=speeds_long |> mutate(dhp=as.POSIXct(str_sub(day_h,4,20),format="%HH%M")) |>
  mutate(dhp=lubridate::ceiling_date(dhp, "15 minutes")) |>
  group_by(link_id,day,dhp) |> 
  summarise(speed=mean(speed,na.rm=TRUE),ts=first(day_h))
speed_temp = speed_quarter |> filter(day=="D9",hour(dhp)>=7,dhp<=as.POSIXct("10H30",format="%HH%M"))
ts.list = unique(speed_temp$ts)
ttf=length(ts.list)
nb = sf::st_relate(st_as_sf(shenzhen.net.pick,"edges") ,st_as_sf(shenzhen.net.pick,"edges") , pattern = "F***T****")
class(nb)="list"

###############################
### Clustering temp 
# indep + match
plots=list()
plots_dendos=list()
cl_list=list()
cut_height=5
for (cts in ts.list){
  
  v = speed_temp |> filter(ts==cts) |> select(speed)

  shenzhen_temp=bind_cols(links.sf,v)
  
  shenzhen.net.temp = as_sfnetwork(shenzhen_temp)
  shenzhen.net.temp = convert(shenzhen.net.temp, to_spatial_smooth,summarise_attributes=list(speed=\(col){mean(col,na.rm=TRUE)}))
  
  Xs = st_as_sf(shenzhen.net.temp,"edges")
  if(sum(is.na(Xs$speed))>0){
    Xs$speed[is.na(Xs$speed)]=sapply(nb[is.na(Xs$speed)],\(nei){mean(Xs$speed[nei],na.rm=TRUE)})
  }
  hc=gtclust_lines(Xs |>select(speed) |>mutate(speed=3.6*speed),method = gtmethod_bayes_dgmm(beta = 20,kappa=1,tau=0.0001))
  plots_dendos[[cts]]=plot(hc)
  cl_list[[cts]]=cutree(hc,h=cut_height)
  
  gg= Xs |> bind_cols(tibble(cl=factor(cutree(hc,h = cut_height))))
  tv=gg |> group_by(cl) |> summarise(vs=var(speed,na.rm = TRUE),n=n()) |> mutate(tv=vs*n) |> filter(n>2) |> summarise(tv=sum(tv)) /(nrow(gg)*var(gg$speed))
  plots[[cts]] = ggplot(gg)+geom_sf(aes(color=cl),show.legend = FALSE)+theme_void()+ggtitle(paste(cts,"TV:",round(tv$tv*100)/100,"K:",length(levels(gg$cl))))
}

ggarrange(plotlist = plots,ncol=5,nrow=3)
ggarrange(plotlist = plots_dendos,ncol=5,nrow=3)

library(forcats)
# matching based on length of intersection
clust_temp = data.frame(matrix(unlist(cl_list),length(cl_list[[1]]),ttf))
colnames(clust_temp)=unique(speed_temp$ts)
nbc = max(clust_temp[,1])
for (icts in 2:length(ts.list)){
  print(icts)
  cl_prec=factor(clust_temp[,(icts-1)])
  cl_cur=factor(clust_temp[,icts])
  
  net_prec = st_as_sf(shenzhen.net.temp,"edges") |> 
    bind_cols(cl_prec=cl_prec) |>
    group_by(cl_prec) |>
    summarise()
  net_cur = st_as_sf(shenzhen.net.temp,"edges") |> 
    bind_cols(cl_cur=cl_cur)|>
    group_by(cl_cur) |>
    summarise()
  
  Le_mat=st_intersection(net_prec,net_cur) |> 
    mutate(L=st_length(geometry)) |> 
    select(cl_prec,cl_cur,L) |> 
    st_drop_geometry() |> 
    tidyr::pivot_wider(names_from=cl_cur,values_from = L,values_fill = 0)
  row_max=apply(Le_mat[,-1],2,which.max)
  col_max=apply(Le_mat[,-1],1,which.max)
  exist_in_previous=col_max[row_max]==(1:(ncol(Le_mat)-1))
  nb_new=sum(!exist_in_previous)
  iperm = as.character(1:(ncol(Le_mat)-1)) 
  iperm[exist_in_previous]=as.character(Le_mat$cl_prec[row_max[exist_in_previous]])
  iperm[!exist_in_previous]=as.character((nbc+1):(nbc+nb_new))
  
  iperm
  nbc = nbc+nb_new
  new_levels=levels(cl_cur)
  names(new_levels)=iperm
  print(new_levels)
  clust_temp[,icts]=fct_recode(cl_cur,!!!new_levels)
}

clust_temp = clust_temp |> mutate_all(\(col){factor(col,levels=as.character(1:nbc))})





# matching based on link count
# clust_temp = data.frame(matrix(unlist(cl_list),nbr,ttf))
# colnames(clust_temp)=unique(speed_temp$ts)
# nbc = max(clust_temp[,1])
# for (icts in 2:length(ts.list)){
#   print(icts)
#   cl_prec=factor(clust_temp[,(icts-1)])
#   cl_cur=factor(clust_temp[,icts])
#   cl_table=table(cl_prec,cl_cur)
#   print(cl_table)
#   iperm = levels(cl_prec)[apply(cl_table,2,which.max)]
#   iperm
#   iprob_prec = apply(cl_table,2,max)/table(cl_prec)[iperm]
#   iprob_prec
#   iprob_cur = apply(cl_table,2,max)/table(cl_cur)
#   iprob_cur
#   iprob=iprob_prec
#   iperm[iprob<0.5]=as.character((nbc+1):(nbc+sum(iprob<0.5)))
#   iperm
#   nbc = nbc+sum(iprob<0.5)
#   new_levels=levels(cl_cur)
#   names(new_levels)=iperm
#   print(new_levels)
#   clust_temp[,icts]=fct_recode(cl_cur,!!!new_levels)
#   print(table(cl_prec,clust_temp[,icts]))
# }
# 
# clust_temp = clust_temp |> mutate_all(\(col){factor(col,levels=as.character(1:nbc))})
plots=list()

d3cat20=c("#ff7f0e","#1f77b4","#2ca02c","#d62728","#9467bd","#ffbb78","#8c564b","#aec7e8","#98df8a","#ff9896","#c5b0d5","#c49c94","#dbdb8d","#bcbd22","#e377c2","#f7b6d2","#7f7f7f","#c7c7c7","#17becf","#9edae5")

for (cts in ts.list){
  v = speed_temp |> filter(ts==cts) |> select(speed)
  
  shenzhen_temp=bind_cols(links.sf,v)
  
  shenzhen.net.temp = as_sfnetwork(shenzhen_temp)
  shenzhen.net.temp = convert(shenzhen.net.temp, to_spatial_smooth,summarise_attributes=list(speed=\(col){mean(col,na.rm=TRUE)}))
  
  Xs = st_as_sf(shenzhen.net.temp,"edges")
  if(sum(is.na(Xs$speed))>0){
    Xs$speed[is.na(Xs$speed)]=sapply(nb[is.na(Xs$speed)],\(nei){mean(Xs$speed[nei],na.rm=TRUE)})
  }
  gg=Xs |> bind_cols(tibble(cl=clust_temp[[cts]]))
  tv=gg |> group_by(cl) |> summarise(vs=var(speed),n=n()) |> mutate(tv=vs*n) |> filter(n!=1) |> summarise(tv=sum(tv)) /(nrow(gg)*var(gg$speed))
  
  str_time = format(speed_temp |> filter(ts==cts) |> pull(dhp) |> unique(),"%H:%M")
  plots[[cts]] = ggplot(gg)+geom_sf(aes(color=cl),show.legend = FALSE)+
    scale_color_manual(drop=FALSE,values=d3cat20)+
    theme_void()+
    ggtitle(paste(str_time,", TV:",round(tv$tv*100)/100))
}
ggarrange(plotlist = plots,ncol=5,nrow=3)



##############################"
## Clustering / temp / graphe => network + time => no mean change + ! topology

ttf=length(ts.list)
nbr = length(nb)
nb_temp.list=lapply(0:(ttf-1),\(tf){lapply(1:length(nb),\(i){
  cnei = nb[[i]]+tf*nbr
  if(tf>0){
    cnei=c(cnei,(tf-1)*nbr+i)
  }
  if(tf<(ttf-1)){
    cnei=c(cnei,(tf+1)*nbr+i)
  }
})})
nb_temp = do.call(c,nb_temp.list)
length(nb_temp)
speeds.list=list()
for(tf in unique(speed_temp$ts)){
  shenzhen_temp=bind_cols(links.sf,speed_temp |> filter(ts==tf) |>  ungroup()|>select(speed))
  shenzhen.net.temp = as_sfnetwork(shenzhen_temp)
  shenzhen.net.temp = convert(shenzhen.net.temp, to_spatial_smooth,summarise_attributes=list(speed=\(s){mean(s,na.rm=TRUE)}))
  speeds.list[[tf]]=st_as_sf(shenzhen.net.temp,"edges") |> st_drop_geometry() |> select(speed)
}
Xs=do.call(rbind,speeds.list)
nrow(Xs)
sum(is.na(Xs))
Xs$speed[is.na(Xs)]=sapply(nb_temp[is.na(Xs)],\(nei){mean(Xs$speed[nei],na.rm=TRUE)})
sum(is.na(Xs))


hc_temp=gtclust_graph(nb_temp,Xs|>mutate(speed=3.6*speed),method = "ward",scaling="raw")


clust_temp = data.frame(matrix(cutree(hc_temp,12),nbr,ttf)) |> mutate_all(\(col){factor(col,levels=1:12)})
colnames(clust_temp)=unique(speed_temp$ts)

shenzhen.cl = bind_cols(st_as_sf(shenzhen.net.temp,"edges"),clust_temp)

sm <- matrix(c(2, 1.2, 0, 1), 2, 2)

shenzhen.cl_tilt <- shenzhen.cl
shenzhen.cl_tilt$geometry <- shenzhen.cl$geometry * sm

library(ggplot2)
library(patchwork)

mpal="Paired"
plots <- lapply(unique(speed_temp$ts)[1:5], function(vname) {
  p <- ggplot(shenzhen.cl) +
    geom_sf(aes_string(color = vname), show.legend = FALSE) +
    scale_color_brewer(palette=mpal)+
    ggtitle(vname)+
    theme_void()
  p
})

p_r1=purrr::reduce(plots, `|`)

plots2 <- lapply(unique(speed_temp$ts)[6:10], function(vname) {
  p <- ggplot(shenzhen.cl) +
    geom_sf(aes_string(color = vname), show.legend = FALSE) +
    scale_color_brewer(palette=mpal)+
    ggtitle(vname)+
    theme_void()
  p
})


p_r2=purrr::reduce(plots2, `|`)


plots3 <- lapply(unique(speed_temp$ts)[11:15], function(vname) {
  p <- ggplot(shenzhen.cl) +
    geom_sf(aes_string(color = vname), show.legend = FALSE) +
    scale_color_brewer(palette=mpal)+
    ggtitle(vname)+
    theme_void()
  p
})

p_r3=purrr::reduce(plots3, `|`)


p_r1 / p_r2 / p_r3
