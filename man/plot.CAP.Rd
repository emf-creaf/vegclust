\encoding{UTF-8}
\name{plot.CAP}
\alias{plot.CAP}
\alias{plot.stratifiedvegdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draws cummulative abundance profiles
}
\description{
Create plots used to inspect one or more cumulative abundance profiles.
}
\usage{
\method{plot}{CAP}(x, sizes=NULL, species=NULL, plots=NULL, switchAxes=FALSE, 
                   add=FALSE, drawAxes = TRUE, xlab="", ylab="", type="s",...)
\method{plot}{stratifiedvegdata}(x, sizes=NULL, species=NULL, plots=NULL, switchAxes=FALSE, 
                   add=FALSE, drawAxes = TRUE, xlab="", ylab="", type="s",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object returned from function \code{\link{CAP}} or an object of class \code{stratifiedvegdata} (see documentation for function \code{\link{stratifyvegdata}}).}
  \item{sizes}{A vector containing the size values associated to each size class. If \code{NULL} the y-axis will be defined using the size class order in \code{x}.}
  \item{species}{A vector of strings indicating the species whose profile is to be drawn. If \code{NULL} all species are plotted.}
  \item{plots}{A vector indicating the plot records whose profile is to be drawn. Can be a \code{character} vector (for plot names), a \code{numeric} vector (for plot indices) or a \code{logical} vector (for TRUE/FALSE selection). If \code{NULL} all plot records are plotted.}
  \item{switchAxes}{A flag indicating whether ordinate and abscissa axes should be interchanged.}
  \item{add}{A flag indicating whether profiles should be drawn on top of current drawing area. If \code{add=FALSE} a new plot is created.}
  \item{drawAxes}{A flag indicating whether axes should be drawn.}
  \item{xlab}{String label for the x axis.}
  \item{ylab}{String label for the y axis.}
  \item{type}{Type of plot to be drawn ("p" for points, "l" for lines, "s" for steps, ...).}
  \item{...}{Additional plotting parameters.}
}
\references{
De \enc{Cáceres}{Caceres}, M., Legendre, P. & He, F. (2013) Dissimilarity measurements and the size structure of ecological communities. Methods in Ecology and Evolution 4: 1167-1177.

}
\author{
Miquel De \enc{Cáceres}{Caceres}, CREAF
}
\seealso{
\code{\link{CAP}}
}
\examples{
## Load stratified data
data(medreg)

## Check that 'medreg' has correct class
class(medreg)

## Create cumulative abundance profile (CAP) for each plot
medreg.CAP = CAP(medreg)

## Draw the stratified data and profile corresponding to the third plot
plot(medreg, plots="3")
plot(medreg.CAP, plots="3")

## Look at the plot and CAP of the same plot
medreg[["3"]]
medreg.CAP[["3"]]
}
