library(R.matlab)
library(dplyr)
library(sf)
library(ggplot2)
library(gtclust)

data("shenzen")

ggplot(shenzen)+geom_sf(aes(color=speed*3.6),size=1.1)+
  scale_color_distiller("Speed (km/h)",palette="RdYlGn",direction = 1,limits=c(0,65))+
  theme_void()+
  theme(legend.position="left")+
  ggtitle("Shenzen speed distribution","8h30-8h45")

cols =c("#c57000",
                  "#003aad",
                  "#74ba00",
                  "#bf51e7",
                  "#929d00",
                  "#718bff",
                  "#8adb78",
                  "#ff75d5",
                  "#008b4f",
                  "#810013",
                  "#b1abff",
                  "#e4c561",
                  "#2f1f58",
                  "#ff9085",
                  "#4a1049",
                  "#cd94ce")


shenzen.df = shenzen |> 
  select(speed)

hc_res=gtclust_lines(shenzen.df,gtmethod_bayes_dgmm())
plot(hc_res)


themebx=theme(panel.grid.minor = element_blank(),
              axis.ticks.y = element_blank(),        ## <- this line
              axis.text.y = element_blank(),         ## <- and this line
              panel.grid.major = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank())

gt_res= geocutree(hc_res,3)
fig1 = ggplot(gt_res)+
  geom_sf(aes(color=factor(cl)))+
  theme_void()+
  scale_color_manual(guide="none",values=cols)

fig2=ggplot(gt_res)+
  geom_sf(aes(color=speed*3.6))+
  theme_void()+
  scale_color_distiller("Speed (km/h)",palette="RdYlGn",direction = 1,limits=c(0,70),guide="none")+
  theme(legend.position="left")+
  ggtitle("K=3")
ggarrange(fig1,fig2)

gg_c=shenzen.df |> bind_cols(cl=cutree(hc_res,3)) 
figbox1  = ggplot(gg_c)+
  geom_boxplot(aes(x=speed*3.6,y=factor(cl),color=factor(cl)))+
  theme_bw()+themebx+
  scale_color_manual(guide="none",values=cols)+labs(x="Speed (Km/h)",y="")

gt_res= geocutree(hc_res,9)
fig3 = 
  ggplot(gt_res)+
  geom_sf(aes(color=factor(cl)))+
  theme_void()+
  scale_color_manual(guide="none",values=cols)


fig4=ggplot(gt_res)+
  geom_sf(aes(color=speed*3.6))+
  theme_void()+
  scale_color_distiller("Speed (km/h)",palette="RdYlGn",direction = 1,limits=c(0,70),guide="none")+
  theme(legend.position="left")+
  ggtitle("K=9")
ggarrange(fig3,fig4)

gg_c=shenzen.df |> bind_cols(cl=cutree(hc_res,9)) 
figbox2  = ggplot(gg_c)+
  geom_boxplot(aes(x=speed*3.6,y=factor(cl),color=factor(cl)))+
  theme_bw()+themebx+
  scale_color_manual(guide="none",values=cols)+labs(x="Speed (Km/h)",y="")

gt_res= geocutree(hc_res,16)
fig5 = ggplot(gt_res)+
  geom_sf(aes(color=factor(cl)))+
  theme_void()+
  scale_color_manual(guide="none",values=cols)

fig6=ggplot(gt_res)+
  geom_sf(aes(color=speed*3.6))+
  theme_void()+
  scale_color_distiller("Speed (km/h)",palette="RdYlGn",direction = 1,limits=c(0,70),guide="none")+
  theme(legend.position="left")+
  ggtitle("K=16")
ggarrange(fig5,fig6)

gg_c=shenzen.df |> bind_cols(cl=cutree(hc_res,16)) 
figbox3  = ggplot(gg_c)+
  geom_boxplot(aes(x=speed*3.6,y=factor(cl),color=factor(cl)))+
  theme_bw()+themebx+
  scale_color_manual(guide="none",values=cols)+labs(x="Speed (Km/h)",y="")

figbox3
ggarrange(plotlist = list(fig6,fig4,fig2,fig5,fig3,fig1,figbox3,figbox2,figbox1))






gg75=shenzen.df |> bind_cols(cl=cutree(hc_res,9)) 
gt_res75 = geocutree(hc_res,9) |> mutate(L=st_length(geometry)) |> arrange(desc(L)) 
mainclusts = gt_res75 |> filter(L>2000) |> pull(cl)

cl_levs = c(-1,mainclusts)

gg75 |> mutate(cl = factor(if_else(!cl %in% mainclusts,-1,cl),levels=cl_levs)) |> ggplot() + geom_boxplot(aes(x=speed*3.6,y=cl,color=cl))+scale_color_discrete(guide="none")+theme_bw()
gt_res75 = gt_res75 |> mutate(cl = if_else(!cl %in% mainclusts,-1,cl))
ggplot(gt_res75)+
  geom_sf(aes(color=factor(cl,cl_levs)))+
  theme_void()+
  theme(legend.position="left")+
  ggtitle("Shenzen clustering results","8h30-8h45")

clust_buff = lapply(1:9,\(ck){
  concaveman::concaveman(st_centroid(gg75|> filter(cl==ck))) |> st_buffer(50)
})

clust_boundaries = st_sf(code_clust=1:9,geom=st_as_sfc(do.call(rbind,clust_buff),crs=NA))

ggplot(shenzen)+
  geom_sf(data=clust_boundaries,color="black",alpha=0)+
  geom_sf(aes(color=speed*3.6),size=1.1)+
  scale_color_distiller("Speed (km/h)",palette="RdYlGn",direction = 1,limits=c(0,65))+
  theme_void()+
  theme(legend.position="left")+
  ggtitle("Shenzen speed distribution","8h30-8h45")
