\encoding{UTF-8}
\name{vegclust-package}
\alias{vegclust-package}
\docType{package}
\title{
Functions for fuzzy and hard clustering of vegetation data
}
\description{
This package contains functions used to perform fuzzy and hard clustering of vegetation data under different models.
}
\details{
\tabular{ll}{
Package: \tab vegclust\cr
Type: \tab Package\cr
Version: \tab 1.6.4\cr
Date: \tab 2016-05-19\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

This package is intended to be helpful for generating and maintaining classifications of vegetation. There are two main functions: (1) \code{\link{vegclust}} is used to create a fuzzy classification from an initial data set (a set of releves of community composition); (2) \code{\link{vegclass}} is used to classify new data (i.e. new releves) into previous classification structures.
}
\author{
Miquel De \enc{Cáceres}{Caceres}, Forest Science Center of Catalonia

Maintainer: Miquel De \enc{Cáceres}{Caceres} <miquelcaceres@gmail.com>
}
\references{
Forgy, E. W. (1965) Cluster analysis of multivariate data: efficiency vs interpretability of classifications. Biometrics 21, 768-769.

MacQueen, J. (1967) Some methods for classification and analysis of multivariate observations. In Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability, eds L. M. Le Cam and J. Neyman, 1, pp. 281-297. Berkeley, CA: University of California Press.

\enc{Davé}{Dave}, R. N. and R. Krishnapuram (1997) Robust clustering methods: a unified view. IEEE Transactions on Fuzzy Systems 5, 270-293.

Bezdek, J. C. (1981) Pattern recognition with fuzzy objective functions. Plenum Press, New York.

Krishnapuram, R. and J. M. Keller. (1993) A possibilistic approach to clustering. IEEE transactions on fuzzy systems 1, 98-110.

De \enc{Cáceres}{Caceres}, M., Font, X, Oliva, F. (2010) The management of numerical vegetation classifications with fuzzy clustering methods. Journal of Vegetation Science 21 (6): 1138-1151.
}
\keyword{ package }
\examples{
## Loads data  
data(wetland)
  
## This equals the chord transformation 
## (see also  \code{\link{decostand}} in package 'vegan')
wetland.chord = as.data.frame(sweep(as.matrix(wetland), 1, 
                              sqrt(rowSums(as.matrix(wetland)^2)), "/"))

## Create noise clustering with 3 clusters. Perform 10 starts from random seeds 
## and keep the best solution
wetland.nc = vegclust(wetland.chord, mobileCenters=3, m = 1.2, dnoise=0.75, 
                      method="NC", nstart=10)

}
