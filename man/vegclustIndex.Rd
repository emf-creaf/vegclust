\encoding{UTF-8}
\name{vegclustIndex}
\alias{vegclustIndex}
\title{
Compute fuzzy evaluation statistics
}
\description{
Computes several evaluation statistics on the fuzzy clustering results on objects of class \code{\link{vegclust}}. 
}
\usage{
vegclustIndex(y)
}
\arguments{
  \item{y}{
	An object of class \code{\link{vegclust}} or a membership matrix.
  }
}
\details{
These statistics were conceived to be computed on fuzzy partitions, such as the ones coming from Fuzzy C-means (Bezdek 1981). Maximum values of PCN or minimum values of PEN can be used as criteria to choose the number of clusters. 
}
\value{
Returns an vector of four values: partition coefficient (PC), normalized partition coefficient (PCN), partition entropy (PE) and normalized partition entropy (PEN).
}
\references{
Bezdek, J. C. (1981) Pattern recognition with fuzzy objective functions. Plenum Press, New York.
}
\author{
Miquel De \enc{Cáceres}{Caceres}, Forest Science Center of Catalonia
}

\seealso{
\code{\link[e1071]{cmeans}},\code{\link{vegclust}}
}
\examples{

## Loads data  
data(wetland)
  
## This equals the chord transformation 
## (see also  \code{\link{decostand}} in package vegan)
wetland.chord = as.data.frame(sweep(as.matrix(wetland), 1, 
                              sqrt(rowSums(as.matrix(wetland)^2)), "/"))
                              
## Create noise clustering with 2, 3 and 4 clusters. Perform 10 starts from random seeds 
## and keep the best solutions
wetland.fcm2 = vegclust(wetland.chord, mobileCenters=2, m = 1.2, method="FCM", nstart=10)
wetland.fcm3 = vegclust(wetland.chord, mobileCenters=3, m = 1.2, method="FCM", nstart=10)
wetland.fcm4 = vegclust(wetland.chord, mobileCenters=4, m = 1.2, method="FCM", nstart=10)

## Compute statistics. Both PCN and PEN indicate that three groups are more advisable 
## than 2 or 4.
print(vegclustIndex(wetland.fcm2))
print(vegclustIndex(wetland.fcm3))
print(vegclustIndex(wetland.fcm4))

}
