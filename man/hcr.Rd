\encoding{UTF-8}
\name{hcr}
\alias{hcr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Heterogeneity-constrained random resampling (HCR) 
}
\description{
Returns a set of indices of the original data set that maximizes the mean and minimizes the variance of the distances between pairs of plot records.
}
\usage{
hcr(d, nout, nsampl=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
	An object of class \code{\link{dist}} containing the distance values between pairs of sites (plot records).
   }
  \item{nout}{The number of sites (plot records) to be chosen among those available in \code{d}.}
  \item{nsampl}{The number of resampling trials to be compared.}
}
\details{
Many subsets of the input data are selected randomly. These subsets are sorted by decreasing mean dissimilarity between pairs of plot records, and then sorted again by increasing variance of these dissimilarities. Ranks from both sortings are summed for each subset, and the subset with the lowest summed rank is considered as the most representative.
}
\value{
Returns a vector containing the indices of the selected sites (plot records) to be used for sub-setting the original table.
}
\references{
Lengyel, A., Chytry, M., Tichy, L. (2011) Heterogeneity-constrained random resampling of phytosociological databases. Journal of Vegetation Science 22: 175-183.
}
\author{
Miquel De \enc{Cáceres}{Caceres}, CREAF
}
\examples{
## Loads data (38 columns and 33 species)
data(wetland)
dim(wetland)
  
## Constructs the chord distance matrix
## (see also  \code{\link{decostand}} in package vegan)
wetland.chord =dist(as.data.frame(sweep(as.matrix(wetland), 1, 
                    sqrt(rowSums(as.matrix(wetland)^2)), "/")))

## Performs HCR resampling. Returns indices of objects
sel = hcr(wetland.chord, nout=20, nsampl=1000)

## Prints the names of the plot records
print(row.names(wetland)[sel])

## Subset the original distance matrix
sel.chord = as.dist(as.matrix(wetland.chord)[sel,sel])
}
