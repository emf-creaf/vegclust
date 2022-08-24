\encoding{UTF-8}
\name{incr.vegclust}
\alias{incr.vegclust}
\title{
Noise clustering with increasing number of clusters
}
\description{
Performs several runs of function 'vegclust' on a community data matrix using an increasing number of clusters until some conditions are met.
}
\usage{
incr.vegclust(x, method="NC", ini.fixed.centers = NULL, 
              min.size = 10, max.var=NULL, alpha = 0.5, 
              nstart=100, fix.previous = TRUE, dnoise=0.75, m=1.0,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	Community data table. A site (rows) by species (columns) matrix or data frame.
  }
  \item{method}{
  A clustering model. Current accepted models are of the noise clustering family: 
  \itemize{
    \item{\code{"NC"}:}{ Noise clustering (Dave and Krishnapuram 1997)}
    \item{\code{"NCdd"}:}{ Noise clustering with medoids}
    \item{\code{"HNC"}:}{ Hard noise clustering}
    \item{\code{"HNCdd"}:}{ Hard noise clustering with medoids}
  }  
  }
  \item{ini.fixed.centers}{
	The coordinates of initial fixed cluster centers. These will be used as \code{fixedCenters} in all calls to \code{\link{vegclust}}. If \code{method="NCdd"} or \code{method="HNCdd"} then \code{ini.fixed.centers} can be specified as a vector of indices for medoids.
  }
  \item{min.size}{
	 The minimum size (cardinality) of clusters. If any of the current k clusters does not have enough members the algorithm will stop and return the solution with k-1 clusters.
   }
  \item{max.var}{
   The maximum variance allowed for clusters (see function \code{\link{clustvar}}). If specified, the algorithm will stop when any of the clusters is at the same time small and has large variance. If \code{max.var = NULL} then this criterion is not used.
   }
  \item{alpha}{
	Criterion to choose cluster seeds from the noise class. Specifically, an object is considered as cluster seed if the membership to the noise class is larger than \code{alpha}.
   }
  \item{nstart}{
	A number indicating how many random trials should be performed for number of groups. Each random trial uses the k-1 cluster centers plus the coordinates of the current cluster seed as initial solution for \code{\link{vegclust}}. Thus, if there are less cluster seed candidates than \code{nstart}, then not all runs are conducted.
  }
  \item{fix.previous}{
 	Flag used to indicate that the cluster centers found when determining k-1 clusters are fixed when determining k clusters.
  }
  \item{m}{
  The fuzziness exponent.
  }
  \item{dnoise}{
	The distance to the noise cluster. 
  }
  \item{...}{
	Additional parameters for function \code{\link{vegclust}}.
  }
}
\details{
Function \code{hier.vegclust} takes starting cluster configurations from cuts of a dendrogram given by object \code{hclust}. Function \code{random.vegclust} chooses random objects as cluster centroids and for each number of clusters performs \code{nstart} trials.
}
\value{
Returns an object of class \code{\link{vegclust}}; or \code{NULL} if the initial cluster does not contain enough members.
}
\author{
Miquel De \enc{Cáceres}{Caceres}, CREAF
}

\references{
\enc{Davé}{Dave}, R. N. and R. Krishnapuram (1997) Robust clustering methods: a unified view. IEEE Transactions on Fuzzy Systems 5, 270-293.
}
\seealso{
\code{\link{vegclust}},\code{\link{hier.vegclust}}
}
\examples{
## Loads data  
data(wetland)
  
## This equals the chord transformation 
## (see also \code{\link{decostand}} in package 'vegan')
wetland.chord = as.data.frame(sweep(as.matrix(wetland), 1,
                              sqrt(rowSums(as.matrix(wetland)^2)), "/"))

## Call incremental noise clustering 
wetland.nc = incr.vegclust(wetland.chord, method="NC", m = 1.2, dnoise=0.75, 
                           min.size=5)

## Inspect cluster sizes
print(wetland.nc$size)

}
