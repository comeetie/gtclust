
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GTclust : A package for fast clustering of spatial or temporal data with contiguity constrained hierarchical clustering

<!-- badges: start -->
<!-- badges: end -->

GTclust builds on top of `?gtclust_graph` it’s main function to offers
fast clustering of spatial or temporal data with contiguity constrained
hierarchical clustering. `?gtclust_graph` is a quite classical
hierarchical clustering but which is designed to takes advantage of
contiguity constraints defined with a graph between data-points. The
contiguity naturally create a sparsely connected graph that can be
leveraged to speed-up the calculations, thanks to efficient
data-structure (Ambroise et al. 2019), from
![O(N^2D)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;O%28N%5E2D%29 "O(N^2D)")
to
![O(M(log(M)+D))](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;O%28M%28log%28M%29%2BD%29%29 "O(M(log(M)+D))"),
with
![M](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;M "M")
the number of links in the contiguity graph. This allow to process 10^5
data-points in seconds. To reach these performances, the dissimilarity
matrix can’t be computed and the algorithm resort to a more simple
“stored data” approach (Murtagh and Contreras 2012), which even if known
to be less efficient in the general case are well fitted for the
contiguity constrained problems. This approach allows to have a low
spatial complexity of
![O(N)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;O%28N%29 "O(N)")
but is less efficient with `single` and `complete` linkage criterions.
Furthermore, this make the use of non-euclidean dissimilarity measures
impossible. Still one may used the classical linkage criterion
compatible with a storage based approach and available in `hclust` or
`agnes` :

-   `ward` minimum within-cluster variance
-   `centroid` or WPGMC
-   `median` or UPGMC

Furthermore, GTclust also offers two Bayesian linkage criterion (Heller
and Ghahramani 2005) that enable model selection :

-   `bayes_mom` mixture of Multinomial for counts data
-   `bayes_dgmm` diagonal Gaussian mixture models for continuous
    features

To ease, the contiguity graph creation process, gtclust offers several
interfaces to easily works with geographical, temporal (the gt in
GTclust comes from here) or sequential data:

-   `?gtclust_temp` to cluster sequential data, the contiguity graph
    follow from the data ordering
-   `?gtclust_poly` to cluster data associated to geographical polygons,
    the contiguity graph follow from shared boundaries
-   `?gtclust_delaunay` to cluster data associated to geographical
    points, the contiguity graph is derived from the Delaunay
    triangulation of the points (thanks to the `? RTriangle` package
    (Shewchuk 1996))
-   `?gtclust_knn` to cluster data associated to geographical points,
    the contiguity graph is derived from the symmetrized knn graph of
    the geographical points (thanks to the `?RANN` package (Arya et al.
    2019))
-   `?gtclust_dist` to cluster data associated to geographical points,
    the contiguity graph is derived from a threshold over distance the
    geographical points

It also offers several methods dedicated ease the manipulations of
results with spatial data thanks to a tight integration with the sf
package (Pebesma 2018).

## Installation

You can install the development version of gtclust from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("comeetie/gtclust")
```

## Example

This is a basic example, we first prepare some spatial polygons data,
here we used modes share in several municipalities around Paris.

``` r
library(gtclust)
library(dplyr)
library(sf)
library(ggplot2)
library(ggpubr)
data("modesshare.idf")

# compute row sums
modesshare.idf <- modesshare.idf |> 
  rowwise(CODE_IRIS) |> 
  mutate(total = sum(c_across(nodep:tcom)))

# normalize per rows
modesshare.idf.percent = modesshare.idf |> 
  filter(total!=0) |>
  transmute(across(nodep:tcom,\(v){v/total})) 
  
```

To do the clustering, we use the `poly` flavor of GTclust and just
provide the polygons data.frame we had just prepared. Then we may use
the classical function from `?hclust`, `?cutree` to cut the dendrogram
at a specific level of the `?plot.gtclust` method to draw the
dendrogram:

``` r
hc=gtclust_poly(modesshare.idf.percent,method="ward")
#> Warning: Some features were not numeric and have been removed from the
#> clustering.
plot(hc)+
  ggtitle("Dendrogram of the modesshare.idf dataset","with ward linkage")+
  scale_y_continuous("Within-cluster intertia")
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
cutree(hc,k=30) |> head(20)
#>  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 
#>  1  1  2  3  4  5  6  4  7  3  8  9  1  1 10 11  4  8 12  4
```

<img src="man/figures/README-clustering-1.png" width="100%" />

You may also use the `?geocutree` function which build directly a
spatial data.frame with the clustering results:

``` r
modesshare_agg = geocutree(hc,k=500)
```

<img src="man/figures/README-plot-1.png" width="100%" />

    #> Warning: attribute variables are assumed to be spatially constant throughout all
    #> geometries

<img src="man/figures/README-plotparis-1.png" width="100%" />

``` r
res_mom = gtclust_delaunay(modesshare.idf |> st_centroid(),method="bayes_mom")
ggplot(data.frame(ts=res_mom$test.stat,K=(nrow(res_mom$data)-1):1))+
    geom_point(aes(x=K,y=ts))+theme_bw()+scale_y_continuous("Test statistic")+ggtitle("Bayesian test score","with respect to the number of cluster")
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

<img src="man/figures/README-benchplot-1.png" width="100%" />

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Ambroise2019" class="csl-entry">

Ambroise, Christophe, Alia Dehman, Pierre Neuvial, Guillem Rigaill, and
Nathalie Vialaneix. 2019. “Adjacency-Constrained Hierarchical Clustering
of a Band Similarity Matrix with Application to Genomics.” *Algorithms
for Molecular Biology* 14 (1): 22.

</div>

<div id="ref-RANN" class="csl-entry">

Arya, Sunil, David Mount, Samuel E. Kemp, and Gregory Jefferis. 2019.
*RANN: Fast Nearest Neighbour Search (Wraps ANN Library) Using L2
Metric*. <https://CRAN.R-project.org/package=RANN>.

</div>

<div id="ref-Heller2005" class="csl-entry">

Heller, Katherine A., and Zoubin Ghahramani. 2005. “Bayesian
Hierarchical Clustering.” In *Proceedings of the 22nd International
Conference on Machine Learning*, 297–304. ICML ’05. New York, NY, USA:
Association for Computing Machinery.
<https://doi.org/10.1145/1102351.1102389>.

</div>

<div id="ref-Murtagh2012" class="csl-entry">

Murtagh, Fionn, and Pedro Contreras. 2012. “Algorithms for Hierarchical
Clustering: An Overview.” *WIREs Data Mining and Knowledge Discovery* 2
(1): 86–97. https://doi.org/<https://doi.org/10.1002/widm.53>.

</div>

<div id="ref-Pebesma2018" class="csl-entry">

Pebesma, Edzer. 2018. “<span class="nocase">Simple Features for R:
Standardized Support for Spatial Vector Data</span>.” *The R Journal* 10
(1): 439–46. <https://doi.org/10.32614/RJ-2018-009>.

</div>

<div id="ref-RTriangle" class="csl-entry">

Shewchuk, Jonathan Richard. 1996. “Triangle: Engineering a 2d Quality
Mesh Generator and Delaunay Triangulator.” In *Applied Computational
Geometry: Towards Geometric Engineering*, edited by Ming C. Lin and
Dinesh Manocha, 1148:203–22. Lecture Notes in Computer Science.
Springer-Verlag.

</div>

</div>
