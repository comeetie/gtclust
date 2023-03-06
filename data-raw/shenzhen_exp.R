
library(dplyr)
library(sfnetworks)
library(sf)
library(ggplot2)
library(ggpubr)
library(gtclust)

data("shenzhen.net")

my_pal = "Set3"
col_back = '#aaaaaa'
col_back_ct = '#555555'
shenzen = st_as_sf(shenzhen.net,"edges")

map_speed = ggplot(shenzen)+geom_sf(aes(color=speed*3.6),size=0.6)+
  scale_color_distiller("Speed (km/h)",palette="RdYlGn",direction = 1,limits=c(10,60),guide=guide_colourbar(ticks.colour = "black"))+
  theme_void()+theme(legend.position=c(0.2,0.8),panel.background = element_rect(fill = col_back, color = col_back_ct))
map_speed

shenzen.df = shenzen |> 
  select(speed) 


hc_res=gtclust_lines(shenzen.df,gtmethod_bayes_dgmm())
fig_dendo=plot(hc_res)+ggtitle("")
fig_dendo
Ka = hc_res$Kunif
Ka=5
gt_res= geocutree(hc_res,Ka)

theme_blank=theme(panel.grid.minor = element_blank(),
                  axis.ticks.y = element_blank(),        
                  axis.text.y = element_blank(),         
                  panel.grid.major = element_blank(),
                  panel.border = element_blank(),panel.background = element_blank())


map_clust = ggplot(gt_res)+
  geom_sf(aes(color=factor(cl)),size=0.6)+
  theme_void()+theme(panel.background = element_rect(fill = col_back, color = col_back_ct))+
  scale_color_brewer(guide="none",palette=my_pal)

gg_c= shenzen.df |> bind_cols(cl=cutree(hc_res,Ka)) 

fig_box  = ggplot(gg_c)+
  geom_boxplot(aes(x=speed*3.6,y=factor(cl),fill=factor(cl)))+
  theme_bw()+theme_blank+
  scale_fill_brewer(guide="none",palette=my_pal)+labs(x="Speed (Km/h)",y="")

map_clust

ggarrange(plotlist = list(map_speed,map_clust,fig_dendo,fig_box),heights =c(0.7,0.3))
ggsave("./data-raw/images/shenzhen.pdf",width=8.5,height=8.5)

