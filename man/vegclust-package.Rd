\encoding{UTF-8}
\name{vegclust-package}
\alias{vegclust-package}
\alias{vegclust}
\docType{package}
\title{
\packageTitle{vegclust}
Functions for fuzzy and hard clustering of vegetation data
}
\description{
\packageDescription{vegclust}
This package contains functions used to perform fuzzy and hard clustering of vegetation data under different models.
}
\details{
The DESCRIPTION file:
\packageDESCRIPTION{vegclust}
\packageIndices{vegclust}
}
\author{
\packageAuthor{vegclust}
 Maintainer: \packageMaintainer{vegclust} 
}
\references{
De Caceres, M., Font, X, Oliva, F. (2010) The management of numerical vegetation classifications with fuzzy clustering methods. Journal of Vegetation Science 21 (6): 1138-1151.

De Cáceres, M., Legendre, P., & He, F. 2013. Dissimilarity measurements and the size structure of ecological communities (D. Faith, Ed.). Methods in Ecology and Evolution 4: 1167–1177.
}
\keyword{package}
\examples{
## Loads data  
data(wetland)
  
## This equals the chord transformation 
wetland.chord = as.data.frame(sweep(as.matrix(wetland), 1, 
                              sqrt(rowSums(as.matrix(wetland)^2)), "/"))

## Create noise clustering with 3 clusters. Perform 10 starts from random seeds 
## and keep the best solution
wetland.nc = vegclust(wetland.chord, mobileCenters=3, m = 1.2, dnoise=0.75, 
                      method="NC", nstart=10)
}