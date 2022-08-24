\encoding{UTF-8}
\name{interclustdist}
\alias{interclustdist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the distance between pairs of cluster centroids 
}
\description{
Calculates the distance between pairs of cluster centroids, given a distance matrix and a cluster vector.
}
\usage{
interclustdist(x, cluster)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	A site-by-site data matrix or an object of class \code{\link{dist}} containing the distance values between pairs of sites (plot records).
   }
  \item{cluster}{
	A vector indicating the hard membership of each object in \code{x} to a set of groups. Can contain \code{NA} values.
   }
}
\value{
An object of class \code{\link{dist}} containing the distances between cluster centers.
}
\author{
Miquel De \enc{Cáceres}{Caceres}, CREAF
}
\examples{
##TO BE DONE##
}
