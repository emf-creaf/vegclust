\encoding{UTF-8}
\name{CAP}
\alias{CAP}
\alias{CAP2matrix}
\alias{CAPcenters}
\alias{CAPquantile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cumulative abundance profile (CAP)
}
\description{
Functions to calculate cumulative abundance profiles (CAPs), to build matrices from them, and to summarize several profiles.
}
\usage{
CAP(x, transform=NULL, verbose=FALSE)
CAP2matrix(CAP, type="cumulative",classWeights=NULL)
CAPcenters(CAP, y=NULL)
CAPquantile(CAP, q = 0.5, y = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A stratified vegetation data set (see function \code{\link{stratifyvegdata}}).}
  \item{transform}{A function or the name of a function to be applied to each cumulative abundance value.}
  \item{verbose}{A logical flag to indicate extra output.}
  \item{CAP}{An object of class '\code{CAP}'.}
  \item{type}{The type of information that the resulting matrix should contain. Either \code{"profile"}, \code{"abundance"} or \code{"volume"}.}
  \item{classWeights}{A numerical vector containing the weight for size class. If \code{NULL}, then all classes are assumed to have the same weight.}
  \item{y}{A vector used as a factor to calculate average or quantile profiles per each level. Alternatively, an object of class \code{\link{vegclust}} for which CAP centroids or medoids are desired.}
  \item{q}{Probability value for which the quantile is desired. By default the median is given.}
}
\details{
Function \code{CAP} replaces the abundance value of a size class by the sum of abundances in this and larger size classes (strata). Thus, upper classes contain smaller abundance values than lower classes, creating a cumulative abundance profile. Function \code{CAP2matrix} takes an object of class '\code{CAP}' and returns a data matrix, where values differ depending on parameter \code{type}: (1) \code{type="cumulative"} simply reshapes the '\code{CAP}' object (a list) into a matrix with as many rows as plot records and where columns are organized in blocks (there are as many blocks as species and each block has as many columns as size classes); (2) \code{type="total"} returns a plot-by-species matrix where each value is the total abundance of the species in the plot (i.e. the CAP value at the ground level); (3) \code{type="volume"} returns a plot-by-species matrix where each value is the sum of CAP values across size classes (a measure of the "volume" occupied by the species in the plot). When provided, \code{classWeights} are used to weight size classes of the cumulative abundance profiles (for (1) and (3) only). Function \code{CAPcenters} calculates the average abundance profile for a set of plot records. If \code{y} is a factor, it is used to speficy groups of samples for which average profiles are to be calculated. If \code{y} is an object of class '\code{\link{vegclust}}' then the function returns the CAP centroids or medoids corresponding to the clustering result. Function \code{CAPquantile} calculates a quantile profile for a set of CAPs. The usage of \code{y} is the same as for \code{CAPcenters}.
}
\value{
Function \code{CAP} returns an object of class '\code{CAP}', similar to objects of class '\code{stratifiedvegdata}' but where abundance values of upper size classes have beed added to those of lower size classes. Function \code{CAP2matrix} returns a matrix with species as rows (columns depend on the value of \code{type}). Functions \code{CAPcenters} and \code{CAPquantile} return an object of class '\code{CAP}'.
}
\references{
De \enc{Cáceres}{Caceres}, M., Legendre, P. & He, F. (2013) Dissimilarity measurements and the size structure of ecological communities. Methods in Ecology and Evolution 4: 1167-1177.

De \enc{Cáceres}{Caceres}, M., Coll, L., \enc{Martín-Alcón}{Martin-Alcon}, S., \enc{González-Olabarria}{Gonzalez-Olabarria}, J.R. (submitted) A general method for the classification of forest stands using structure and composition.
}
\author{
Miquel De \enc{Cáceres}{Caceres}, CREAF.
}
\seealso{
  \code{\link{stratifyvegdata}}, \code{\link{plot.CAP}}, \code{\link{vegdiststruct}}
}
\examples{
## Load stratified data
data(medreg)

## Check that 'medreg' has correct class
class(medreg)

## Look at the data for the third plot
medreg[[3]]

## Create cumulative abundance profile (CAP) for each plot
medreg.CAP = CAP(medreg)

## Look at the profile of the third plot
medreg.CAP[[3]]

## Create matrix with species abundances
medreg.X = CAP2matrix(medreg.CAP, type="total")
head(medreg.X)

## Generate and plot average profile
average.CAP = CAPcenters(medreg.CAP)
plot(average.CAP)

## Generate and plot median profile
median.CAP = CAPquantile(medreg.CAP, q = 0.5)
plot(median.CAP)

}
