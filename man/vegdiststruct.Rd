\encoding{UTF-8}
\name{vegdiststruct}
\alias{vegdiststruct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Structural and compositional dissimilarity
}
\description{
Function to calculate the dissimilarity between ecological communities taking into account both their composition and the size of organisms.
}
\usage{
vegdiststruct(x, y=NULL, paired=FALSE, type="cumulative", method="bray", 
              transform=NULL, classWeights=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A stratified vegetation data set (see function \code{\link{stratifyvegdata}}), a set of cummulative abundance profiles (see function \code{\link{CAP}}) or a set of cummulative abundance surfaces (see function \code{\link{CAS}}).}
  \item{y}{A second stratified vegetation data set (see function \code{\link{stratifyvegdata}}), a second set of cummulative abundance profiles (see function \code{\link{CAP}}) or a second set of cummulative abundance surfaces (see function \code{\link{CAS}}) against which object \code{x} should be compared.}
  \item{paired}{Only relevant when \code{y != NULL}. If \code{paired = TRUE} pairwise comparisons are calculated between elements in \code{x} and \code{y} (and \code{x} and \code{y} need to be of the same length). If \code{paired = FALSE} then all objects in \code{x} are compared to all objects in \code{y}.}
  \item{type}{Whether dissimilarities between pairs of sites should be calculated from differences in cummulative abundance (\code{"cumulative"}), in total abundance (\code{"total"}) or in volumes of cumulative abundance profiles (\code{"volume"}).}
  \item{method}{The dissimilarity coefficient to calculate (see details).}
  \item{transform}{A function or the name of a function to be applied to each cumulative abundance value.}
  \item{classWeights}{A numerical vector or a matrix containing the weight of each size class or combination of size classes (see functions \code{\link{CAP2matrix}} and \code{\link{CAS2matrix}}). If \code{NULL}, then the function assumes classes of equal weight.}
}
\details{
The six different coefficients available are described in De Caceres et al. (2013): (1) \code{method="bray"} for percentage difference (alias Bray-Curtis dissimilarity); (2) \code{method="ruzicka"} for Ruzicka index (a generalization of Jaccard); (3) \code{method="kulczynski"} for the Kulczynski dissimilarity index; (4) \code{method="ochiai"} for the complement of a quantitative generalization of Ochiai index of similarity; (5) \code{method="canberra"} for the Canberra index (Adkins form); (6) \code{method="relman"} for the relativized Manhattan coefficient (Whittaker's index of association). Currently, the function also supports (7) \code{method="manhattan"} for the city block metric.
}
\value{
Returns an object of class '\code{\link{dist}}'.
}
\references{
De \enc{Cáceres}{Caceres}, M., Legendre, P. & He, F. (2013) Dissimilarity measurements and the size structure of ecological communities. Methods in Ecology and Evolution 4: 1167-1177.
}
\author{
Miquel De \enc{Cáceres}{Caceres}, Forest Science Center of Catalonia.
}
\seealso{
\code{\link{stratifyvegdata}}, \code{\link[vegan]{vegdist}}
}
\examples{
## Load stratified data
data(medreg)

## Check that 'medreg' has correct class
class(medreg)

## Create cumulative abundance profile (CAP) for each plot
medreg.CAP = CAP(medreg)

## Create dissimilarity (percentage difference) matrix using profiles
medreg.D = vegdiststruct(medreg, method="bray")

## Create dissimilarity (percentage difference) matrix using abundances
medreg.D2 = vegdiststruct(medreg, method="bray", type="total")

## Calculate correlation
cor(as.vector(medreg.D), as.vector(medreg.D2))

}
