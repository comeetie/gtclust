library(gtools)
library(ggplot2)
n=250
#pi=0.6+rnorm(n)*0.01
#lpi=c(sum(log(pi)),sum(log(1-pi)))

P = rdirichlet(n,c(35,200))
d=ncol(P)
lpi = colSums(log(P))
lambda=c(0.1,0.1)
Cn=function(alpha){lgamma(sum(alpha))-sum(lgamma(alpha))}
pdf_r = function(beta){
  sum(exp(beta)*(lpi-lambda))+sum(log(lambda)+beta)+n*Cn(exp(beta)+1)
}




# derivative
# ! lpi = sum(log(pi))
g_reppostdiririchlet=function(beta,lpi,n,lambda){
  alpha=exp(beta)
  alpha*(n*digamma(sum(alpha+1))-n*digamma(alpha+1)+(lpi-lambda))+1
}


gd_reppostdiririchlet=function(beta,lpi,n,lambda){
  -exp(2*beta)*n*trigamma(exp(beta)+1)
}

z_reppostdiririchlet=function(beta,lpi,n,lambda){
  n*trigamma(sum(exp(beta)+1))
}


up_reppostdiririchlet=function(beta,lpi,n,lambda){
  g  = g_reppostdiririchlet(beta,lpi,n,lambda)
  gd = gd_reppostdiririchlet(beta,lpi,n,lambda)
  z  = z_reppostdiririchlet(beta,lpi,n,lambda)
  Di = diag(1/(g+gd-1))
  nf = (1/z+sum(exp(2*beta)/(g+gd-1)))
  off_diag = outer(exp(beta),exp(beta))
  off_diag = off_diag/(g+gd-1)
  off_diag = t(t(off_diag)/(g+gd-1))
  Hi = Di-1/nf*off_diag
  as.vector(Hi%*%g)
}




ldet_reppostdiririchlet=function(beta,lpi,n,lambda){
  g  = g_reppostdiririchlet(beta,lpi,n,lambda)
  gd = gd_reppostdiririchlet(beta,lpi,n,lambda)
  z  = z_reppostdiririchlet(beta,lpi,n,lambda)
  log(1+z*sum(exp(2*beta)/(g+gd-1)))+sum(log(g+gd-1))
}

reppostdiririchlet=function(beta0,lpi,n,lambda,alpha=0.5,bd = 0.5){
  beta = beta0
  nbit = 5000
  i=0
  delta=10
  while(delta>10^-8 & i<nbit){
    i=i+1
    beta_old=beta
    g = g_reppostdiririchlet(beta,lpi,n,lambda)
    v = up_reppostdiririchlet(beta,lpi,n,lambda)
    t = 1
    beta_temp = beta-t*v
    # while(pdf_r(beta_temp)<(pdf_r(beta)+alpha*t*sum(g*v))){
    #   t=bd*t
    #   beta_temp = beta-t*v
    # }
    beta = beta_temp
    delta = pdf_r(beta)-pdf_r(beta_old)
    print(delta)
  }
  print(i)
  beta
}



xy_opt = optim(c(0,0),\(x){-pdf_r(x)})

bb1=c(xy_opt$par[1]-2.5,xy_opt$par[1]+2.5)
bb2=c(xy_opt$par[2]-2.5,xy_opt$par[2]+2.5)

xy=expand.grid(seq(bb1[1],bb1[2],length.out=500),seq(bb2[1],bb2[2],length.out=500))
lmpi=log(colMeans(P))
pdf=function(alpha){apply(alpha,1,pdf_r)}
xyz=cbind(xy,pdf(xy))
pdf.surf=data.frame(alpha1=xyz[,1],alpha2=xyz[,2],pdf=xyz[,3])
ggplot(pdf.surf)+
  geom_contour_filled(aes(x=alpha1,y=alpha2,z=pdf))+
  geom_point(data=data.frame(x=xy_opt$par[1],y=xy_opt$par[2]),aes(x=x,y=y),color="red")+
  geom_abline(aes(slope=1,intercept=xy_opt$par[2]-xy_opt$par[1]),color="red")+
  geom_abline(aes(slope=1,intercept=lmpi[2]-lmpi[1]),color="green")+
  scale_color_distiller(palette = "Spectral", direction = -1)



s0 = n*(d)

alpha0 = s0*colMeans(P)

beta0 = log(alpha0-1)

beta_hat = reppostdiririchlet(beta0,lpi,n,lambda)

exp(beta_hat)+1
xy_opt$par

ggplot(pdf.surf)+
  geom_contour_filled(aes(x=alpha1,y=alpha2,z=exp(pdf)))+
  geom_point(data=data.frame(x=xy_opt$par[1],y=xy_opt$par[2]),aes(x=x,y=y),color="red")+
  geom_abline(aes(slope=1,intercept=xy_opt$par[2]-xy_opt$par[1]),color="red")+
  geom_point(data=data.frame(x=lmpi[1]+log(10),y=lmpi[2]+log(10)),aes(x=x,y=y),color="green")+
  geom_abline(aes(slope=1,intercept=lmpi[2]-lmpi[1]),color="green")+
  geom_point(data=data.frame(x=beta0[1],y=beta0[2]),aes(x=x,y=y),color="blue")+
  geom_abline(aes(slope=1,intercept=beta0[2]-beta0[1]),color="blue")+
  geom_point(data=data.frame(x=beta_hat[1],y=beta_hat[2]),aes(x=x,y=y),color="black")+
  scale_color_distiller(palette = "Spectral", direction = -1)

