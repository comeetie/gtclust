N = 100000
K = 50
D = 3
clusters_sizes = rpois(K-1,N/K)
i_change = c(0,cumsum(clusters_sizes),N)
X = matrix(0,nrow=N,ncol=D)
cl = rep(0,N)
for (k in 2:(K+1)) {
  ind = (i_change[k-1]+1):i_change[k]
  nk = length(ind)
  cl[ind]=k-1
  #X[ind,]=cbind(rnorm(nk,runif(1)*5),rnorm(nk,runif(1)*5),rnorm(nk,runif(1)*5))
  X[ind,]=cbind(rpois(nk,runif(1)*10),rpois(nk,runif(1)*10),rpois(nk,runif(1)*10))
}


sol=gtclust_temp(X,method="bayes_mom")
clh = cutree(sol,50)
image(table(cl,clh))
