
benchmark_ades <- function(nobj) {

  res <- matrix(NA,length(nobj),5) |> as.data.frame()
  colnames(res) <- c("N.objects","Storage","Time","M.links","Dist.time")
  res[,1L] <- nobj

  for(i in 1:length(nobj)) { 
    N <- nobj[i]
    coords.mem <- cbind(x=runif(N,-1,1),y=runif(N,-1,1))
    dat.mem <- runif(N,0,1)
    if(i>1L) rm(D.mem) ; gc()
    {start.time = Sys.time()
      D.mem <- try(dat.mem %>% dist)
    end.time = Sys.time()}
    res[i,5L] <- end.time-start.time
    if(any(class(D.mem)=="try-error"))
      break
    neighbors.mem <-
      (coords.mem %>%
         tri2nb %>%
         nb2listw(style="B") %>%
         listw2sn)[,1:2]
    {start.time = Sys.time()
      res.mem <- try(constr.hclust(D.mem, method="ward.D2",
                                   neighbors.mem))
      end.time = Sys.time()}
    if(any(class(res.mem)=="try-error"))
      break
    res[i,2L] <- (2*object_size(D.mem) + object_size(neighbors.mem) +
                    object_size(res.mem))/1048576  # n. bytes per MiB
    res[i,3L] <- end.time-start.time
    res[i,4L] <- nrow(neighbors.mem)
    
  }
  res[["N.objects"]] <- as.integer(res[["N.objects"]])
  res[["M.links"]] <- as.integer(res[["M.links"]])
  res
}

res_adespatial <- benchmark_ades(nobj=c(1000,2000,5000,10000,20000,50000,100000))


benchmark_geoh <- function(nobj) {
  # Argument -
  # nobj : Number of objects in simulation runs
  res <- matrix(NA,length(nobj),4) |> as.data.frame()
  colnames(res) <- c("N.objects","Storage","Time","M.links")
  res[,1L] <- nobj
  ## resources <- list()
  for(i in 1:length(nobj)) { 
    N <- nobj[i]
    coords.mem <- cbind(x=runif(N,-1,1),y=runif(N,-1,1))
    dat.mem <- matrix(runif(N,0,1),ncol=1)
    neighbors.mem <- tri2nb(coords.mem)
       
    {start.time = Sys.time()
      res.mem <- try(geohclust_graph(neighbors.mem,dat.mem))
      end.time = Sys.time()}
    if(any(class(res.mem)=="try-error"))
      break
    res[i,2L] <- (object.size(neighbors.mem) +
                    object.size(res.mem))/1048576  # n. bytes per MiB
    res[i,3L] <- end.time-start.time
    res[i,4L] <- sum(card(neighbors.mem))
  }
  res[["N.objects"]] <- as.integer(res[["N.objects"]])
  res[["M.links"]] <- as.integer(res[["M.links"]])
  res
}
res_geoh <- benchmark_geoh(nobj=c(1000,2000,5000,10000,20000,50000,100000))
res_geoh[["Dist.time"]]=0

library(ggplot2)
ggplot()+geom_line(data=res_geoh,aes(x=N.objects,y=Time),col="red")+
  geom_point(data=res_geoh,aes(x=N.objects,y=Time),col="red")+
  geom_line(data=res_adespatial,aes(x=N.objects,y=Dist.time),col="blue")+
  geom_point(data=res_adespatial,aes(x=N.objects,y=Dist.time),col="blue")

ggplot()+geom_line(data=res_geoh,aes(x=N.objects,y=Storage),col="red")+
  geom_point(data=res_geoh,aes(x=N.objects,y=Storage),col="red")+
  geom_line(data=res_adespatial,aes(x=N.objects,y=Storage),col="blue")+
  geom_point(data=res_adespatial,aes(x=N.objects,y=Storage),col="blue")

