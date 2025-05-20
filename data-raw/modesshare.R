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

download.file("https://wxs.ign.fr/1yhlj2ehpqf3q6dt6a2y7b64/telechargement/inspire/CONTOURS-IRIS-2020-01-01$CONTOURS-IRIS_2-1__SHP__FRA_2020-01-01/file/CONTOURS-IRIS_2-1__SHP__FRA_2020-01-01.7z","./data-raw/contours-iris.7z",method="wget")
archive_extract(archive= "./data-raw/contours-iris.7z",dir = "./data-raw/")

iris=read_sf("./data-raw/CONTOURS-IRIS_2-1__SHP__FRA_2020-01-01/CONTOURS-IRIS/1_DONNEES_LIVRAISON_2020-12-00282/CONTOURS-IRIS_2-1_SHP_LAMB93_FXX-2020/CONTOURS-IRIS.shp")

modesshare.idf.raw = iris |> left_join(data.act,by=c("CODE_IRIS"="IRIS")) |>
  filter(!is.na(nodep)) |>
  mutate(DEP=substr(CODE_IRIS,1,2)) |> 
  filter(DEP %in% c(75,77,78,91,92,93,95,94)) 

library(rmapshaper)

modesshare.idf = ms_simplify(modesshare.idf.raw,  keep=0.6) |> st_make_valid() |> st_cast("MULTIPOLYGON")

usethis::use_data(modesshare.idf,overwrite = TRUE)

modesshare.rc.raw = iris |> left_join(data.act,by=c("CODE_IRIS"="IRIS")) |>
  filter(!is.na(nodep)) |>
  mutate(DEP=substr(CODE_IRIS,1,2)) |> 
  filter(DEP %in% c(18,28,36,37,41,45))

modesshare.rc = ms_simplify(modesshare.rc.raw,  keep=0.6) |> st_make_valid() |> st_cast("MULTIPOLYGON")

usethis::use_data(modesshare.rc,overwrite = TRUE)

modesshare.pts = iris |> st_centroid() |> left_join(data.act,by=c("CODE_IRIS"="IRIS")) 

usethis::use_data(modesshare.pts,overwrite = TRUE)

system("rm -rf ./data-raw/CONTOURS-IRIS_2-1__SHP__FRA_2021-01-01")
system("rm -rf ./data-raw/CONTOURS-IRIS_2-1__SHP__FRA_2020-01-01")
system("rm -rf ./data-raw/contours-iris.7z")
system("rm -rf ./data-raw/act.zip")
system("rm -rf ./data-raw/base-ic-activite-residents-2018.CSV")
system("rm -rf ./data-raw/meta_base-ic-activite-residents-2018.CSV")
system("rm -rf ./data-raw/meta_base-ic-activite-residents-2018.CSV")

