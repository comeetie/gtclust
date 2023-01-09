library(sf)
library(dplyr)
library(ggplot2)
library(gtclust)
library(mclust)
library(rgeoda)

square=st_polygon(list(rbind(c(0,0),c(0,100),c(100,100),c(100,0),c(0,0))))
psize = 30
K = 9
Kr = sqrt(K)
clust_size = psize/Kr

grid=st_make_grid(square,n=c(psize,psize))
nb=st_intersects(grid,grid)
class(nb)="list"
cl=do.call(c,lapply(1:Kr,\(r){
rep(rep(c(((r-1)*Kr+1):(r*Kr)),each=clust_size),clust_size)
}))


grid.sf=st_as_sf(grid)
grid.sf$cl=factor(cl)

ggplot(grid.sf)+
  geom_sf(aes(fill=cl),color="#000000")+
  scale_fill_brewer(palette="Set3",guide="none")+
  theme_void()



centers = c(1,5,2,3,9,7,8,6,4)
matrix(centers,3,3)

sigma=1.8
grid.sf = grid.sf |> mutate(v=rnorm(nrow(grid.sf),centers[cl],sigma))

qw=rook_weights(grid.sf)




sigma_list=c(0.25,0.5,0.75,1,1.25,1.5)
res=list()
for (sigma in sigma_list){
  clist = list()
  for (r in 1:50){
    grid.sf=grid.sf |> mutate(v=rnorm(nrow(grid.sf),centers[cl],sigma))
    sol=gtclust_poly(grid.sf|>select(v),adjacency = "queen",method = gtclust::gtmethod_bayes_dgmm(),display_progress = TRUE)
    mod=Mclust(grid.sf$v,G = 1:30,modelNames = "V")
    mod9=Mclust(grid.sf$v,G = 9,modelNames = "V")
    df=grid.sf|>st_drop_geometry()|>select(v)
    
    sk=skater(9,qw,df)
    aricode::NMI(grid.sf$cl,sk$Clusters)
    
    azpg = azp_greedy(9, qw, df)
    aricode::NMI(grid.sf$cl,azpg$Clusters)
    
    azpsa = azp_sa(9, qw, df, cooling_rate = 0.85)
    aricode::NMI(grid.sf$cl,azpsa$Clusters)
    redc = redcap(9, qw, df, "fullorder-wardlinkage")
    aricode::NMI(grid.sf$cl,redc$Clusters)
    
    clist[[r]]=rbind(
      data.frame(r=r,sigma=sigma,alg="mclust",nmi=aricode::NMI(grid.sf$cl,mod$classification),K=mod$G),
      data.frame(r=r,sigma=sigma,alg="gtclust",nmi=aricode::NMI(grid.sf$cl,cutree(sol,sol$Kunif)),K=sol$Kunif),
      data.frame(r=r,sigma=sigma,alg="mclust9",nmi=aricode::NMI(grid.sf$cl,mod9$classification),K=9),
      data.frame(r=r,sigma=sigma,alg="gtclust9",nmi=aricode::NMI(grid.sf$cl,cutree(sol,9)),K=9),
      data.frame(r=r,sigma=sigma,alg="skater",nmi=aricode::NMI(grid.sf$cl,sk$Clusters),K=9),
      data.frame(r=r,sigma=sigma,alg="azpg",nmi=aricode::NMI(grid.sf$cl,azpg$Clusters),K=9),
      data.frame(r=r,sigma=sigma,alg="azpsa",nmi=aricode::NMI(grid.sf$cl,azpsa$Clusters),K=9),
      data.frame(r=r,sigma=sigma,alg="redcap",nmi=aricode::NMI(grid.sf$cl,redc$Clusters),K=9))
  }
  res[[paste0("sigma_",sigma)]]=clist
}


res.df=do.call(rbind,lapply(res,\(lr){do.call(rbind,lr)}))
res.df.mean = res.df |> group_by(sigma,alg) |>
  summarize(nmi_mean=mean(nmi),K_mean=mean(K))


ggplot(res.df)+
  geom_boxplot(aes(x=alg,y=nmi))+
  facet_wrap(~sigma)+
  coord_flip()+
  theme_bw()+labs(title="NMI with simulated partition",subtitle = "for several values of sigma",y="NMI",x="")

# library(ggtext)
#  labs(title="Mutual information with simulated partition",
#       subtitle="for <span style='color:#66c2a5'>gtclust</span> and <span style='color:#fc8d62'>mclust</span>, on the ten region datasets.",
#       x="sigma",y="NMI")+
#  theme(plot.subtitle = element_markdown())+ylim(c(0,1))








sigma_list=c(0.25,0.75,1.25,1.5)
res=list()
for (sigma in sigma_list){
  clist = list()
  grid.sf=grid.sf |> mutate(v=rnorm(nrow(grid.sf),centers[cl],sigma))
    
  # data
  clist$data=ggplot(grid.sf)+geom_sf(aes(fill=v,color=v))+
      scale_color_distiller(palette="Greys",guide="none")+
      scale_fill_distiller(palette="Greys",guide="none")+
      theme_void()+ggtitle(paste0("Sigma: ",sigma))
    
  # gt
  sol=gtclust_poly(grid.sf|>select(v),adjacency = "rook",method = gtclust::gtmethod_bayes_dgmm(),display_progress = TRUE)
  
  grid.sf$cluster=cutree(sol,sol$Kunif)
  clist$gtclust=ggplot(grid.sf)+
    geom_sf(aes(fill=factor(cluster)),color="#555555",size=0.2)+
    scale_fill_brewer(palette="Set3",guide="none")+
    theme_void()+ggtitle("gtclust")
  
  
  # mclust
  mod=Mclust(grid.sf$v,G = 1:30,modelNames = "V")
  mod9=Mclust(grid.sf$v,G = 9,modelNames = "V")
  grid.sf$cluster=mod$classification
  clist$mclust=ggplot(grid.sf)+
    geom_sf(aes(fill=factor(cluster)),color="#555555",size=0.2)+
    scale_fill_brewer(palette="Set3",guide="none")+
    theme_void()+ggtitle("mclust")
  
  
  df=grid.sf|>st_drop_geometry()|>select(v)
    
  # skater
  sk=skater(9,qw,df)
  grid.sf$cluster=sk$Clusters
  clist$skater=ggplot(grid.sf)+
      geom_sf(aes(fill=factor(cluster)),color="#555555",size=0.2)+
      scale_fill_brewer(palette="Set3",guide="none")+
      theme_void()+ggtitle("skater")
    
    
    
  azpg = azp_greedy(9, qw, df)
  grid.sf$cluster=azpg$Clusters
  clist$azpg=ggplot(grid.sf)+
    geom_sf(aes(fill=factor(cluster)),color="#555555",size=0.2)+
    scale_fill_brewer(palette="Set3",guide="none")+
    theme_void()+ggtitle("azpg")
  res=c(res,clist)
}

ggarrange(plotlist = res)
