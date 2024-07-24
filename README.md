
# Fuzzy clustering of vegetation data <a href="https://emf-creaf.github.io/vegclust/"><img src="man/figures/logo.png" align="right" height="139" alt="vegclust website" /></a>

<!-- badges: start -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/vegclust)](https://cran.r-project.org/package=vegclust)
[![](https://cranlogs.r-pkg.org/badges/vegclust)](https://cran.rstudio.com/web/packages/vegclust/index.html)
[![R-CMD-check](https://github.com/emf-creaf/vegclust/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/emf-creaf/vegclust/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Introduction

Package `vegclust` is a package designed to assist analyses of
vegetation structure and composition. It is intended to be useful for
community ecologists and forest engineers, but the clustering functions
can be used in other fields. The package provides functions to:

- Perform fuzzy clustering of vegetation data (De Cáceres et al. 2010).
- Assess ecological community ressemblance on the basis of structure and
  composition (De Cáceres et al. 2013).

## Package installation

Package vegclust can be found at
[CRAN](https://cran.r-project.org/package=vegclust) package repository.
In addition, the latest stable `vegclust` R package can be installed
from GitHub as follows:

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

## Note about ‘Community Trajectory Analysis’

Until ver. 2.0, package `vegclust` included functions to conduct
Community Trajectory Analysis (CTA). Since ver. 2.0 these functions have
been moved to an independent package `ecotraj` available at
<https://github.com/emf-creaf/ecotraj/> and also in
[CRAN](https://cran.r-project.org/package=ecotraj).

## Maintenance

Although not in active development, the R package is maintained by the
[*Ecosystem Modelling Facility*](https://emf.creaf.cat) unit at CREAF
(in Spain).

## References

- De Cáceres, M., Font, X., & Oliva, F. 2010. The management of
  vegetation classifications with fuzzy clustering. Journal of
  Vegetation Science 21: 1138–1151
  (<https://doi.org/10.1111/j.1654-1103.2010.01211.x>).

- De Cáceres, M., Legendre, P., & He, F. 2013. Dissimilarity
  measurements and the size structure of ecological communities (D.
  Faith, Ed.). Methods in Ecology and Evolution 4: 1167–1177
  (<https://doi.org/10.1111/2041-210X.12116>).
