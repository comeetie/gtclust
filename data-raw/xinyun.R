library(sf)
library(dplyr)
library(jsonlite)
library(gtclust)
library(igraph)
library(patchwork)
data=read_sf('./data-raw/toulouse_sensor_clustering(1).geojson')

data %>% select(detid,flow,k_correct,fclass) %>% st_drop_geometry()

nei = fromJSON("./data-raw/toulouse_sensor_neighbors(1).json")



G = graph_from_adj_list(nei$neighbor,mode="all")
comps = igraph::components(G,mode = "strong")

X = data %>% select(flow,k_correct,fclass) %>% st_drop_geometry() %>% mutate(fclass=factor(fclass))

# composante connexe correspond au 346 premier capteurs
nb=nei$neighbor[1:346,]
X=X[1:346,]
data=data[1:346,]

hc_res = gtclust_graph(nb,X,gtmethod_bayes_mixed())

plot(hc_res)

table(cutree(hc_res,8),X$fclass[1:346])
data$cluster=factor(cutree(hc_res,8))

g1 = ggplot(data) + geom_sf(aes(col=cluster)) + theme_void()
g1

data_agg = data |> st_drop_geometry() |> group_by(cluster) |>
  summarize(mflow=mean(flow),mk=mean(k_correct),vflow=var(flow),vk=var(k_correct),n=n())

data_agg
data_agg_disc = data |> st_drop_geometry() |> group_by(cluster,fclass) |>
  summarize(n=n()) |> 
  group_by(cluster) |> 
  mutate(nt=sum(n)) |> 
  mutate(p=n/nt) |>
  select(cluster,fclass,p) |>
  tidyr::pivot_wider(names_from = fclass,values_from = p,values_fill = 0)

data_agg_disc

g2 = ggplot(data) + 
  geom_point(aes(x=flow,y=k_correct,color=cluster)) + 
  geom_point(data=data_agg,aes(x=mflow,y=mk,color=cluster),size=4,shape = 3)+theme_bw()
g2 

g1|g2

