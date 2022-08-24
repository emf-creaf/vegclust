\encoding{UTF-8}
\name{relate.levels}
\alias{relate.levels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Relates two clustering level results.
}
\description{
Analyzes how lower level clusters are assigned into upper level ones. The analysis is made for several number of clusters.
}
\usage{
relate.levels(lower, upper, defuzzify = FALSE, excludeFixed = FALSE, verbose=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lower}{ A list of objects of type \code{\link{vegclust}} or \code{\link{vegclass}} that represent classifications at a finer level of resolution.}
  \item{upper}{ A list of objects of type \code{\link{vegclust}} or \code{\link{vegclass}} that represent classifications at an broader level of resolution.}
  \item{defuzzify}{ A logical flag used to indicate whether the result of calling \code{\link{crossmemb}} should be deffuzified.}
  \item{excludeFixed}{A logical used to indicate whether fixed clusters should be excluded from the comparison of levels.}
  \item{verbose}{ A flag used to ask for extra screen output.}
  \item{...}{
    Additional parameters for function \code{\link{defuzzify}}.
  }
}
\details{
For each pair of \code{vegclust} (or \code{vegclass}) objects in \code{upper} and \code{lower}, the function calls function \code{\link{crossmemb}} and then, if asked, deffuzifies the resulting memberships (by calling function \code{\link{defuzzify}}) and several quantities are calculated (see 'value' section).
}
\value{
 A list with several data frames (see below). In each of them, the rows are items of \code{upper} and columns are items of \code{lower}. The names of rows and columns are the number of clusters of each \code{\link{vegclust}} (or \code{vegclass}) object.
 \item{nnoise}{The number of low level clusters that are assigned to the Noise class (for \code{upper} objects using Noise clustering). }
 \item{maxnoise}{The maximum membership value of low level clusters to the Noise class (for \code{upper} objects using Noise clustering). }
 \item{minmaxall}{The minimum value (across upper level clusters) of the maximum membership value observed among the lower level clusters. }
 \item{minallsize}{The minimum value (across upper level clusters) of the sum of membership values across lower level clusters. }
 \item{empty}{The number of upper level clusters (mobile or fixed) that do not have any member among the lower level clusters. }
}
\author{
Miquel De \enc{Cáceres}{Caceres}, CREAF
}

\seealso{
\code{\link{vegclust}}, \code{\link{vegclass}}, \code{\link{defuzzify}}
}
\examples{
## Loads data  
data(wetland)
  
## This equals the chord transformation 
## (see also  \code{\link{decostand}} in package vegan)
wetland.chord = as.data.frame(sweep(as.matrix(wetland), 1,
                              sqrt(rowSums(as.matrix(wetland)^2)), "/"))

## Create noise clustering from hierarchical clustering at different number of cluster
wetland.hc = hclust(dist(wetland.chord),method="ward") 
wetland.nc1 = hier.vegclust(wetland.chord, wetland.hc, cmin=2, cmax=6, m = 1.2, 
                            dnoise=0.75, method="NC")
wetland.nc2 = hier.vegclust(wetland.chord, wetland.hc, cmin=2, cmax=4, m = 1.2, 
                            dnoise=0.85, method="NC")

## Studies the assignment of levels
relate.levels(wetland.nc1, wetland.nc2, method="cut")
}
