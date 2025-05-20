library(readr)
library(sf)
library(dplyr)
zones= st_read("../mobitic_fv/data/TRIRIS/Data Transport TRIRIS/TRAG_MOBITIC_zonage_TRIRIS.geojson")|>
  filter(perimetre==1)|> mutate(i=0:(n()-1)) |> st_transform(2154)
ods = read_csv("../mobitic_fv/ods_weekdays.csv") |> 
  mutate(Origine=as.character(Origine),Destination=as.character(Destination)) |>
  inner_join(zones|>st_drop_geometry()|>select(regroupement_desc,io=i),by=c("Origine"="regroupement_desc")) |>
  inner_join(zones|>st_drop_geometry()|>select(regroupement_desc,id=i),by=c("Destination"="regroupement_desc")) |>
  mutate(Volume=round(Volume))

ods

triplet = ods|> select(io,id,Volume) |> as.matrix()


nb = sf::st_relate(zones,zones, pattern = "F***T****")
nb_c = lapply(nb,\(nei){nei-1})




res=gtclust_dcsbm(nb_c,triplet,0.05,0.02)
plot(res$Ll)

# complete inter prior with linear slope if needed
miss_prior = is.na(res$PriorInter)
if(sum(miss_prior)>0){
  nbmiss = max(which(miss_prior))
  res$PriorInter[miss_prior]=seq(res$PriorIntra[length(res$PriorIntra)],res$PriorInter[nbmiss+1],length.out=nbmiss)
}
# compute the spanning tree prior term
ptree = (res$PriorInter+res$PriorIntra-res$PriorInter[1])
Llf = res$Ll + ptree +res$PriorK;
# convert merge mat in hclust format
V = nrow(zones);
merge_mat = apply(res$merge,2,function(col){ifelse(col<V,-(col+1),col-V+1)})
# format the results in hclust form
hc_res = list(merge=merge_mat,
              Ll = -Llf,
              Kunif = length(Llf)-which.max(Llf)+1,
              PriorIntra = res$PriorIntra,
              PriorInter = res$PriorInter,
              PriorK = res$PriorK)


table(cutree(hc_res,k=40))
plot(Llf)


vth=100
ori = ods |>filter(Volume>vth,Origine!="999999",Destination!="999999") |> 
  left_join(zones,by=c("Origine"="regroupement_desc"))|>st_as_sf()|>st_centroid() |> mutate(fid=1:n()) |> 
  select(fid)
dest = ods |>filter(Volume>vth,Origine!="999999",Destination!="999999") |> 
  left_join(zones,by=c("Destination"="regroupement_desc"))|>st_as_sf()|>st_centroid() |> mutate(fid=1:n()) |> 
  select(fid)

links = bind_rows(ori,dest) |> 
  group_by(fid) |> 
  summarise(do_union=FALSE) |> 
  st_cast("LINESTRING") |> 
  mutate(Volume=ods|>filter(Volume>vth,Origine!="999999",Destination!="999999") |>pull(Volume))

ko=nrow(zones)-which.max(Llf)
ko=50
zones_cl = zones |> mutate(cl=as.factor(cutree(hc_res,k=ko))) |> 
  left_join(ods|>count(Origine,wt=Volume,name="VolumeOut"),by=c("regroupement_desc"="Origine")) |>
  left_join(ods|>count(Destination,wt=Volume,name="VolumeIn"),by=c("regroupement_desc"="Destination"))

regions = zones_cl |> group_by(cl) |> summarise(n=n())

plot(regions |> select(cl))

ggplot() + 
  geom_sf(data=regions,aes(fill=factor(cl)),alpha=0.95,linewidth=0)+
  geom_sf(data=links|>filter(Volume>10),aes(linewidth=Volume)) + 
  scale_linewidth_continuous(limits=c(0,200000),range=c(0,8))+
  scale_color_discrete(guide="none")+
  scale_fill_discrete(guide="none")+
  geom_sf(data=regions,color="red",fill="#ffffff00",linewidth=0.4,linetype=2)


ggplot(zones_cl |>st_centroid()) + 
  geom_sf(data=regions,aes(fill=factor(cl)),alpha=0.05,linewidth=0)+
  geom_sf(data=links|>filter(Volume>10),aes(linewidth=Volume)) + 
  scale_linewidth_continuous(limits=c(0,200000),range=c(0,8))+
  geom_sf(aes(color=factor(cl),size=VolumeIn)) + theme_void() + 
  scale_color_discrete(guide="none")+
  scale_fill_discrete(guide="none")+
  geom_sf(data=regions,color="red",fill="#ffffff00",linewidth=0.4,linetype=2)
