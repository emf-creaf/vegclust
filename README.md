vegclust
================

## Introduction

Package `vegclust` is a package designed to assist analyses of
vegetation structure and composition. It is intended to be useful for
community ecologists and forest engineers, but the clustering functions
can be used in other fields. The package provides functions to:

-   Perform fuzzy clustering of vegetation data: De Cáceres et
    al. (2010) (<https://doi.org/10.1111/j.1654-1103.2010.01211.x>).
-   Assess ecological community ressemblance on the basis of structure
    and composition: De Cáceres et al. (2013)
    (<https://doi.org/10.1111/2041-210X.12116>).

## Package installation

Package vegclust can be found at [CRAN](https://cran.r-project.org/)
package repository. In addition, the latest stable `vegclust` R package
can be installed from GitHub as follows:

``` r
devtools::install_github("emf-creaf/vegclust")
```

Additionally, users can have help to run package functions directly as
package vignettes, by forcing their inclusion in installation:

``` r
devtools::install_github("emf-creaf/vegclust", 
                         build_opts = c("--no-resave-data", "--no-manual"),
                         build_vignettes = TRUE)
```

## References

-   De Cáceres, M., Font, X., & Oliva, F. 2010. The management of
    vegetation classifications with fuzzy clustering. Journal of
    Vegetation Science 21: 1138–1151.

-   De Cáceres, M., Legendre, P., & He, F. 2013. Dissimilarity
    measurements and the size structure of ecological communities (D.
    Faith, Ed.). Methods in Ecology and Evolution 4: 1167–1177.
