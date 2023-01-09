centers=c(1,5,7,3,8,2,9,4,6,10)
cl=rep(1:10,each=100)

sigma_list=seq(0.2,4,by=0.2)
res=list()
for (sigma in sigma_list){
  clist = list()
  for (r in 1:10){
    X=data.frame(cl=cl) |> mutate(x=rnorm(length(cl),centers[cl],sigma))
    sol=gtclust_temp(X|>select(x),method = gtclust::gtmethod_bayes_dgmm())
    plot(sol)
    mod=Mclust(X$x,1:20,modelNames = "V")
    clist[[r]]=rbind(
    data.frame(r=r,sigma=sigma,alg="mclust",nmi=aricode::NMI(cl,mod$classification),K=mod$G),
    data.frame(r=r,sigma=sigma,alg="gtclust",nmi=aricode::NMI(cl,cutree(sol,sol$Kunif)),K=sol$Kunif))
  }
  res[[paste0("sigma_",sigma)]]=clist
}


res.df=do.call(rbind,lapply(res,\(lr){do.call(rbind,lr)}))
res.df.mean = res.df |> group_by(sigma,alg) |>
  summarize(nmi_mean=mean(nmi),K_mean=mean(K))

library(ggtext)
ggplot(res.df.mean)+
  geom_point(aes(x=sigma,y=nmi_mean,color=alg))+
  geom_point(data=res.df,aes(x=sigma,y=nmi,color=alg),alpha=0.3)+
  geom_line(aes(x=sigma,y=nmi_mean,color=alg))+
  theme_bw()+
  scale_color_brewer(palette="Set2",guide="none")+
  labs(title="Mutual information with simulated partition",
       subtitle="for <span style='color:#66c2a5'>gtclust</span> and <span style='color:#fc8d62'>mclust</span>, on the ten segment datasets.",
       x="sigma",y="NMI")+
  theme(plot.subtitle = element_markdown())

ggplot(res.df.mean)+
  geom_point(aes(x=sigma,y=K_mean,color=alg))+
  geom_point(data=res.df,aes(x=sigma,y=K,color=alg),alpha=0.3)+
  geom_line(aes(x=sigma,y=K_mean,color=alg))+
  theme_bw()+
  scale_color_brewer(palette="Set2",guide="none")+
  labs(title="Number of clusters",
       subtitle="found by <span style='color:#66c2a5'>gtclust</span> and <span style='color:#fc8d62'>mclust</span>, on the ten segment datasets.",
       x="sigma",y="K")+
  theme(plot.subtitle = element_markdown())
