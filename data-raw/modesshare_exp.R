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

hc=gtclust_poly(modesshare.rc.log |> select(nodep:drm,tcom),method=gtclust::gtmethod_bayes_dgmm())
plot(hc)


ggplot(modesshare.rc)+geom_sf(fill="white",size=0.3)+
  theme_void()


K=23
modesshare.rc.cl = bind_cols(modesshare.rc,cl=cutree(hc,K))
modesshare.rc.cl.small = modesshare.rc.cl |> 
  group_by(cl) |> 
  summarise_if(is.numeric,sum) |>  
  mutate(across(nodep:tcom,\(v){v/total})) 
ggplot(modesshare.rc.cl.small)+geom_sf(fill="white",size=0.3)+
  theme_void()



K=23

labels =c("Tours","Orléans","Bourges","Blois","Chartres","Dreux","Montargis","Vendôme","Vierzon","Châteauroux")

labels.sf=modesshare.rc |> filter(NOM_COM %in% labels) |> 
  group_by(NOM_COM) |> 
  summarise(n=n()) |>
  ungroup() |>
  st_centroid()
library(shadowtext)
xyl=data.frame(labels.sf |> st_coordinates(),label=labels.sf$NOM_COM)

modesshare.rc.cl = bind_cols(modesshare.rc,cl=cutree(hc,K))
modesshare.rc.cl.small = modesshare.rc.cl |> 
  group_by(cl) |> 
  summarise_if(is.numeric,sum) |>  
  mutate(across(nodep:tcom,\(v){v/total})) 


a = ggplot(modesshare.rc.cl.small)+geom_sf(aes(fill=voiture*100),size=0.3)+
  theme_void()+
  geom_shadowtext(data=xyl,aes(x=X,y=Y,label=label),color="#ffffff",fontface = "bold",nudge_y = 8000)+
  scale_fill_distiller("Car (%)",palette="Oranges",direction=1)
a
b = ggplot(modesshare.rc.cl.small)+geom_sf(aes(fill=tcom*100),size=0.3)+
  theme_void()+
  geom_shadowtext(data=xyl,aes(x=X,y=Y,label=label),color="#ffffff",fontface = "bold",nudge_y = 8000)+
  scale_fill_distiller("Transit (%)",palette="Blues",direction=1)
b
c = ggplot(modesshare.rc.cl.small)+geom_sf(aes(fill=marche*100),size=0.3)+
  theme_void()+
  geom_shadowtext(data=xyl,aes(x=X,y=Y,label=label),color="#ffffff",fontface = "bold",nudge_y = 8000)+
  scale_fill_distiller("Walking (%)",palette="Purples",direction=1)
c
d=ggplot(modesshare.rc.cl.small)+geom_sf(aes(fill=(1-(tcom+voiture+marche))*100),size=0.3)+
  theme_void()+
  geom_shadowtext(data=xyl,aes(x=X,y=Y,label=label),color="#ffffff",fontface = "bold",nudge_y = 8000)+
  scale_fill_distiller("Others (%)",palette="Greens",direction=1)
d

ggarrange(b,c)



a2 = ggplot(modesshare.rc.perc)+geom_sf(aes(fill=voiture*100,color=voiture*100),size=0.3)+
  geom_sf(data=modesshare.rc.cl.small,alpha=0)+
  theme_void()+
  geom_shadowtext(data=xyl,aes(x=X,y=Y,label=label),color="#ffffff",fontface = "bold",nudge_y = 8000)+
  scale_fill_distiller("Car (%)",palette="Oranges",direction=1)+
  scale_color_distiller("Car (%)",palette="Oranges",direction=1)
a2

b2 = ggplot(modesshare.rc.log)+geom_sf(aes(fill=tcom,color=tcom),size=0.3)+
  geom_sf(data=modesshare.rc.cl.small,alpha=0)+
  theme_void()+
  geom_shadowtext(data=xyl,aes(x=X,y=Y,label=label),color="#ffffff",fontface = "bold",nudge_y = 8000)+
  scale_fill_distiller("Transit (%)",palette="Blues",direction=1)+
  scale_color_distiller("Transit (%)",palette="Blues",direction=1)
b2

c2 = ggplot(modesshare.rc.log)+geom_sf(aes(fill=marche,color=marche),size=0.3)+
  theme_void()+
  geom_sf(data=modesshare.rc.cl.small,alpha=0)+
  geom_shadowtext(data=xyl,aes(x=X,y=Y,label=label),color="#ffffff",fontface = "bold",nudge_y = 8000)+
  scale_fill_distiller("Walking (%)",palette="Purples",direction=1)+
  scale_color_distiller("Walking (%)",palette="Purples",direction=1)
c2

d2=ggplot(modesshare.rc.perc)+geom_sf(aes(fill=(1-(tcom+voiture+marche))*100,color=(1-(tcom+voiture+marche))*100),size=0.3)+
  theme_void()+
  geom_sf(data=modesshare.rc.cl.small,alpha=0)+
  geom_shadowtext(data=xyl,aes(x=X,y=Y,label=label),color="#ffffff",fontface = "bold",nudge_y = 8000)+
  scale_fill_distiller("Others (%)",palette="Greens",direction=1)+
  scale_color_distiller("Others (%)",palette="Greens",direction=1)
d2


ggbox= modesshare.rc.perc |> bind_cols(cl=cutree(hc,K)) |> st_drop_geometry() |> tidyr::pivot_longer(nodep:tcom,names_to = "mode",values_to = "perc")

ggplot(ggbox)+geom_boxplot(aes(color=mode,y=factor(cl),x=perc))+facet_wrap(~mode)

