% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conformveg.R
\name{conformveg}
\alias{conformveg}
\title{Conform two community data tables}
\usage{
conformveg(x, y, fillvalue = 0, verbose = FALSE)
}
\arguments{
\item{x}{Community data, a site-by-species matrix.}

\item{y}{Community data, a site-by-species matrix.}

\item{fillvalue}{The value to be used to fill new entries in inflated matrices.}

\item{verbose}{Displays information about the number of species shared between \code{x} and \code{y}, as well as the number of species that are in one of the data tables but not in the other.}
}
\value{
A list with the two inflated matrices \code{x} and \code{y}.
}
\description{
Conforms two community data tables to have the same set of columns (species)
}
\details{
This function adds to \code{x} as many new columns as columns of \code{y} that are not in \code{x}. The same is done for \code{y}, so the two tables have the same set of columns when they are returned.
}
\examples{
## Loads data (38 columns and 33 species)
data(wetland)
dim(wetland)

## Splits wetland data into two matrices of 30x27 and 11x22
wetland.30 <- wetland[1:30,]
wetland.30 <- wetland.30[,colSums(wetland.30)>0]
dim(wetland.30)
wetland.11 <- wetland[31:41,]
wetland.11 <- wetland.11[,colSums(wetland.11)>0] 
dim(wetland.11)

## Conforms the two matrices so they can eventually be merged
wetland.cf <- conformveg(wetland.30, wetland.11)
dim(wetland.cf$x)
dim(wetland.cf$y)
names(wetland.cf$x)==names(wetland.cf$y)
}
\seealso{
\code{\link{vegclust}}, \code{\link{vegclass}}
}
\author{
Miquel De \enc{Cáceres}{Caceres}, CREAF.
}
