data(well.df)
library(gtclust)
changepoints = gtclust_temp(well.df,gtmethod_bayes_dgmm(),scaling = "zscore")

plot(changepoints)
nbr=20
ggplot(well.df)+geom_line(aes(x=1:nrow(well.df),y=nmr),color="#555555")+
  geom_vline(aes(xintercept=p),data=data.frame(p=which(diff(cutree(changepoints,nbr))!=0)),color="#ff6666")+
  scale_color_discrete(guide="none")+theme_bw()+labs(x="",title="Change point detection",subtitle = "on the well log dataset")
