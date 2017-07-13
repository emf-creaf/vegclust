\encoding{UTF-8}
\name{vegclust2kmeans}
\alias{vegclust2kmeans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reshapes as kmeans object
}
\description{
This function casts an object of class \code{\link{vegclust}} into an object of class \code{\link{kmeans}}.
}
\usage{
vegclust2kmeans(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{\link{vegclust}} to be casted, where \code{method="KM"} and \code{mode="raw"}.}
}
\author{
Miquel De \enc{Cáceres}{Caceres}, Forest Science Center of Catalonia
}
\seealso{
\code{\link{vegclust}}, \code{\link{kmeans}}
}
\examples{
## Loads data  
data(wetland)
  
## This equals the chord transformation 
## (see also  \code{\link{decostand}} in package vegan)
wetland.chord = as.data.frame(sweep(as.matrix(wetland), 1, 
                              sqrt(rowSums(as.matrix(wetland)^2)), "/"))

## Create noise clustering with 3 clusters. Perform 10 starts from random seeds 
wetland.vc = vegclust(wetland.chord, mobileCenters=3, 
                      method="KM", nstart=10)

## Reshapes as kmeans object
wetland.km = vegclust2kmeans(wetland.vc)
wetland.km
}
