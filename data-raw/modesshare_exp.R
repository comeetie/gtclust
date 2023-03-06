library(dplyr)
library(sf)
library(patchwork)
data("modesshare.rc")

# compute row sums
modesshare.rc <- modesshare.rc |> 
  rowwise(CODE_IRIS) |> 
  mutate(total = sum(c_across(nodep:tcom))) 

imiss=which(modesshare.rc$total<50)
nb = st_intersects(modesshare.rc,modesshare.rc)

for (i in imiss){
  modesshare.rc$nodep[i]=sum(modesshare.rc$nodep[c(nb[[i]],i)])
  modesshare.rc$tcom[i]=sum(modesshare.rc$tcom[c(nb[[i]],i)])
  modesshare.rc$voiture[i]=sum(modesshare.rc$voiture[c(nb[[i]],i)])
  modesshare.rc$drm[i]=sum(modesshare.rc$drm[c(nb[[i]],i)])
  modesshare.rc$marche[i]=sum(modesshare.rc$marche[c(nb[[i]],i)])
  modesshare.rc$velo[i]=sum(modesshare.rc$velo[c(nb[[i]],i)])
  modesshare.rc$total[i]=sum(modesshare.rc$total[c(nb[[i]],i)])
}

# normalize per rows
modesshare.rc.log = modesshare.rc |> 
  group_by(CODE_IRIS,NOM_COM)|>
  transmute(across(nodep:tcom,\(v){log(v+1)-log(voiture+1)})) |> ungroup()

modesshare.rc.perc = modesshare.rc |>
  rowwise(CODE_IRIS) |> 
  mutate(total = sum(c_across(nodep:tcom))) |> 
  group_by(CODE_IRIS,NOM_COM)|>
  transmute(across(nodep:tcom,\(v){v/total})) |> ungroup()

hc=gtclust_poly(modesshare.rc.log |> select(nodep,marche:drm,tcom),method=gtclust::gtmethod_bayes_dgmm(),display_progress = TRUE)
dendo=plot(hc)+ggtitle("")
dendo
ggsave("./data-raw/images/dendo_modesshare.pdf",width=3.5,height=3)

  
K=hc$Kunif

map = geocutree(hc,K)

labels =c("Tours","Orléans","Bourges","Blois","Chartres","Dreux","Montargis","Châteauroux")

labels.sf=modesshare.rc |> filter(NOM_COM %in% labels) |> 
  group_by(NOM_COM) |> 
  summarise(n=n()) |>
  ungroup() |>
  st_centroid()
library(shadowtext)
xyl=data.frame(labels.sf |> st_coordinates(),label=labels.sf$NOM_COM)

modesshare.rc.cl = bind_cols(modesshare.rc,cl=cutree(hc,K))

# mos = modesshare.rc.log |> select(nodep,marche:drm,tcom)
# qw=queen_weights(mos)
# sk=skater(23,qw,mos)
# modesshare.rc.cl = bind_cols(modesshare.rc,cl=sk$Clusters)


modesshare.rc.cl.small = modesshare.rc.cl |> 
  group_by(cl) |> 
  summarise_if(is.numeric,sum) |>  
  mutate(across(nodep:tcom,\(v){v/total})) 


a = ggplot(modesshare.rc.cl.small)+geom_sf(aes(fill=voiture*100),size=0.3)+
  theme_void()+theme(legend.position = "bottom")+
  geom_shadowtext(data=xyl,aes(x=X,y=Y,label=label),color="#ffffff",fontface = "bold",nudge_y = 8000)+
  scale_fill_distiller("Car (%)",palette="Oranges",direction=1)
a
b = ggplot(modesshare.rc.cl.small)+geom_sf(aes(fill=tcom*100),size=0.3)+
  theme_void()+theme(legend.position = "bottom")+
  geom_shadowtext(data=xyl,aes(x=X,y=Y,label=label),color="#dddddd",fontface = "bold",nudge_y = 8000,size=4,family="Helvetica")+
  scale_fill_distiller("Transit (%)",palette="Blues",direction=1)
b
c = ggplot(modesshare.rc.cl.small)+geom_sf(aes(fill=marche*100),size=0.3)+
  theme_void()+theme(legend.position = "bottom")+
  geom_shadowtext(data=xyl,aes(x=X,y=Y,label=label),color="#dddddd",fontface = "bold",nudge_y = 8000,size=4,family="Helvetica")+
  scale_fill_distiller("Walking (%)",palette="Purples",direction=1)
c
d=ggplot(modesshare.rc.cl.small)+geom_sf(aes(fill=(1-(tcom+voiture+marche))*100),size=0.3)+
  theme_void()+theme(legend.position = "bottom")+
  geom_shadowtext(data=xyl,aes(x=X,y=Y,label=label),color="#ffffff",fontface = "bold",nudge_y = 8000)+
  scale_fill_distiller("Others (%)",palette="Greens",direction=1)
d

ggarrange(b,c)
ggsave("data-raw/images/map_modeshare.pdf",width=10,height=6)


modesshare.rc.cl = bind_cols(modesshare.rc,cl=cutree(hc,9))
modesshare.rc.cl.small = modesshare.rc.cl |> 
  group_by(cl) |> 
  summarise_if(is.numeric,sum) |>  
  mutate(across(nodep:tcom,\(v){v/total})) 


a = ggplot(modesshare.rc.cl.small)+geom_sf(aes(fill=voiture*100),size=0.3)+
  theme_void()+theme(legend.position = "bottom")+
  geom_shadowtext(data=xyl,aes(x=X,y=Y,label=label),color="#ffffff",fontface = "bold",nudge_y = 8000)+
  scale_fill_distiller("Car (%)",palette="Oranges",direction=1)
a
b = ggplot(modesshare.rc.cl.small)+geom_sf(aes(fill=tcom*100),size=0.3)+
  theme_void()+theme(legend.position = "bottom")+
  geom_shadowtext(data=xyl,aes(x=X,y=Y,label=label),color="#dddddd",fontface = "bold",nudge_y = 8000,size=4,family="Helvetica")+
  scale_fill_distiller("Transit (%)",palette="Blues",direction=1)
b
c = ggplot(modesshare.rc.cl.small)+geom_sf(aes(fill=marche*100),size=0.3)+
  theme_void()+theme(legend.position = "bottom")+
  geom_shadowtext(data=xyl,aes(x=X,y=Y,label=label),color="#dddddd",fontface = "bold",nudge_y = 8000,size=4,family="Helvetica")+
  scale_fill_distiller("Walking (%)",palette="Purples",direction=1)
c
d=ggplot(modesshare.rc.cl.small)+geom_sf(aes(fill=(1-(tcom+voiture+marche))*100),size=0.3)+
  theme_void()+theme(legend.position = "bottom")+
  geom_shadowtext(data=xyl,aes(x=X,y=Y,label=label),color="#ffffff",fontface = "bold",nudge_y = 8000)+
  scale_fill_distiller("Others (%)",palette="Greens",direction=1)
d

ggarrange(b,c)
ggsave("data-raw/images/map_modeshare_K9.pdf",width=10,height=6)
