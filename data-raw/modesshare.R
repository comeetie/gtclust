#https://www.insee.fr/fr/statistiques/5650708
download.file("https://www.insee.fr/fr/statistiques/fichier/5650708/base-ic-activite-residents-2018_csv.zip","./data-raw/act.zip")
unzip("./data-raw/act.zip",exdir = "./data-raw/")
library(readr)
library(dplyr)
library(archive)
library(sf)
data.act = read_delim("./data-raw/base-ic-activite-residents-2018.CSV",delim = ";") |> 
  select(IRIS,COM,TYP_IRIS,103:108) |>
  mutate(nodep=round(C18_ACTOCC15P_PAS)) |>
  mutate(marche=round(C18_ACTOCC15P_MAR)) |>
  mutate(velo=round(C18_ACTOCC15P_VELO)) |>
  mutate(drm=round(C18_ACTOCC15P_2ROUESMOT)) |>
  mutate(voiture=round(C18_ACTOCC15P_VOIT)) |>
  mutate(tcom=round(C18_ACTOCC15P_TCOM)) |>
  select(1,10:15)

download.file("ftp://Contours_IRIS_ext:ao6Phu5ohJ4jaeji@ftp3.ign.fr/CONTOURS-IRIS_2-1__SHP__FRA_2021-01-01.7z","./data-raw/contours-iris.7z",method="wget")
archive_extract(archive= "./data-raw/contours-iris.7z",dir = "./data-raw/")

iris=read_sf("./data-raw/CONTOURS-IRIS_2-1__SHP__FRA_2021-01-01/CONTOURS-IRIS/1_DONNEES_LIVRAISON_2021-06-00217/CONTOURS-IRIS_2-1_SHP_LAMB93_FXX-2021/CONTOURS-IRIS.shp")

modesshare_idf = iris |> left_join(data.act,by=c("CODE_IRIS"="IRIS")) |>
  filter(!is.na(nodep)) |>
  mutate(DEP=substr(CODE_IRIS,1,2)) |> 
  filter(DEP %in% c(75,77,78,91,92,93,95,94)) 

library(rmapshaper)

modesshare = ms_simplify(modesshare_idf,  keep=0.03) |> st_make_valid() |> st_cast("MULTIPOLYGON")

usethis::use_data(modesshare,overwrite = TRUE)

system("rm -rf ./data-raw/CONTOURS-IRIS_2-1__SHP__FRA_2021-01-01")
system("rm -rf ./data-raw/contours-iris.7z")
system("rm -rf ./data-raw/act.zip")
system("rm -rf ./data-raw/base-ic-activite-residents-2018.CSV")
system("rm -rf ./data-raw/meta_base-ic-activite-residents-2018.CSV")

