\encoding{UTF-8}
\name{crossmemb}
\alias{crossmemb}
\title{
Cross-table of two fuzzy classifications
}
\description{
Calculates a cross-tabulated matrix relating two fuzzy membership matrices
}
\usage{
crossmemb(x, y, relativize = TRUE)
}
\arguments{
  \item{x}{A site-by-group fuzzy membership matrix. Alternatively, an object of class 'vegclust' or 'vegclass'.}
  \item{y}{A site-by-group fuzzy membership matrix. Alternatively, an object of class 'vegclust' or 'vegclass'.}
  \item{relativize}{If \code{TRUE} expresses the cross-tabulated values as proportions of cluster size in \code{x}.}
}
\value{A cross-tabulated matrix comparing the two classifications. In general, each cell's value is the (fuzzy) number of objects that in \code{x} are assigned to the cluster corresponding to the row and in \code{y} are assigned to the cluster corresponding to the column. If \code{relativize=TRUE} then the values of each row are divided by the (fuzzy) size of the corresponding cluster in \code{x}. 
}

\author{
Miquel De \enc{Cáceres}{Caceres}, CREAF.
}

\seealso{
\code{\link{defuzzify}}, \code{\link{vegclust}}
}
\examples{
## Loads data  
data(wetland)
  
## This equals the chord transformation 
## (see also  \code{\link{decostand}} in package vegan)
wetland.chord = as.data.frame(sweep(as.matrix(wetland), 1, 
                              sqrt(rowSums(as.matrix(wetland)^2)), "/"))


## Create clustering with 3 clusters. Perform 10 starts from random seeds 
## and keep the best solution. Try both FCM and NC methods:
wetland.fcm = vegclust(wetland.chord, mobileCenters=3, m = 1.2, method="FCM", nstart=10)
wetland.nc = vegclust(wetland.chord, mobileCenters=3, m = 1.2, dnoise=0.75, method="NC", 
                      nstart=10)

## Compare the results
crossmemb(wetland.fcm, wetland.nc, relativize=FALSE)
}
