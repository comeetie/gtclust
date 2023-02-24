
data(well.df)
library(gtclust)
library(ggplot2)
hc = gtclust_temp(well.df,gtmethod_bayes_dgmm(tau=0.00001,beta = var(well.df$nmr)))

library(latex2exp)
gath =ggplot2::theme(
  panel.grid.major = ggplot2::element_blank(),
  panel.grid.minor = ggplot2::element_blank(),
  panel.border = ggplot2::element_blank(),
  legend.position="bottom",
  plot.title = element_text(family = "Helvetica", face = "bold", size = (15)),
  axis.text = element_text(family = "Helvetica",size=(10)),
  axis.title = element_text(family = "Helvetica", size = (12)))

df=data.frame(Ll=hc$Ll,k=length(hc$Ll):1)
tree=hc
mlogalpha=NULL;
df.line=data.frame(Ll=-tree$Ll[which.min(tree$Ll):length(tree$Ll)],k=(length(tree$Ll)-which.min(tree$Ll)+1):1)
df.front = df.line
if(which.min(tree$Ll)>1){
  df.front$xend = tree$height[(which.min(tree$Ll)-1):length(tree$height)]
}else{
  df.front$xend = c(0,tree$height)
}

df.front = df.front[nrow(df.front):1,]
df.front = df.front[!duplicated(df.front$x),]
if(!is.null(mlogalpha)){
  if(sum(df.front$xend<mlogalpha)>2){
    df.front = df.front[df.front$xend<mlogalpha,]
  }else{
    error("mlogalpha too small, no front to show.",call. = FALSE)
  }
}
xmax = max(df.front$xend)*1.05
ny  = diff(range(df.front$Ll))*0.025
df.front$x=c(xmax,df.front$xend[1:(nrow(df.front)-1)])
df.front$lf = (df.front$x-df.front$xend)/df.front$x[1]

#Fig front

ggplot2::ggplot(df.front)+
  ggplot2::geom_abline(data=df.line,ggplot2::aes_(intercept=~Ll,slope=~-(k-1)),alpha=0.15)+
  ggplot2::geom_abline(data=df.line,ggplot2::aes_(intercept=~Ll,slope=~-(k-1)),color="#cacaca",alpha=0.05)+
  ggplot2::geom_segment(data=df.front,ggplot2::aes_(x=~x,xend=~xend,y=~Ll-x*(k-1),yend=~Ll-xend*(k-1)))+
  ggplot2::geom_point(data=df.front,ggplot2::aes_(x=~xend,y=~Ll-(k-1)*xend),size=1.5)+
  ggplot2::geom_point(data=df.front,ggplot2::aes_(x=~xend,y=~Ll-(k-1)*xend),color="white",size=1)+
  ggplot2::scale_x_continuous(expression(-log(alpha)),expand = expansion(0.01,0))+
  ggplot2::scale_y_continuous(name=TeX("$L^{post}(X,c,K,\\alpha)$"),expand = expansion(0.01,0))+
  ggplot2::theme_bw()+ggtitle("Un-normalized log-posterior",TeX("with respect to $-\\log(\\alpha)$"))+gath
ggsave("./data-raw/images/ex_front.pdf",width=5,height=4)



rownames(tree$data)=1:nrow(tree$data)
if(substr(tree$method,1,5)=="bayes"){
  im = which.min(tree$Ll)-1
  nb_max_leafs=(length(tree$height)+1)-im
}else{
  im=(length(tree$height)+1)-nb_max_leafs
}
small_tree = gtclust:::collapse(tree, nb_max_leafs)
dend_data = ggdendro::dendro_data(as.dendrogram(small_tree))
cluster_sizes = sapply(small_tree$members,length)
dend_data$labels$size=cluster_sizes[small_tree$order]


#Fig dendo

segs_constr = dend_data$segments
ggplot2::ggplot() + 
  ggplot2::geom_segment(data=segs_constr,ggplot2::aes_(y =~ x, x  =~ y, yend =~ xend, xend =~ yend))+
  ggplot2::theme_bw()+
  ggplot2::scale_y_continuous(TeX("$L^{post}(X,c,K,\\alpha)$"),n.breaks = 3,labels = c(10000,20000,30000))+
  geom_vline(data=df.front,aes(xintercept = xend),alpha=0.25)+
  ggplot2::scale_x_continuous(name = expression(-log(alpha)),expand = expansion(0.01,0))+
  theme(
    axis.text.y = element_text(family = "Helvetica",size=(10),color="white"),
    axis.title.y = element_text(family = "Helvetica", size = (12),color="white"),
    axis.ticks.y = element_line(color="white")
  )+gath+
  ggtitle("Dendrogram","with heights derived from log-posterior tipping points.")

ggsave("./data-raw/images/ex_dendo.pdf",width=5,height=4)

