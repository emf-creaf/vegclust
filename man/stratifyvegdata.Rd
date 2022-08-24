\encoding{UTF-8}
\name{stratifyvegdata}
\alias{stratifyvegdata}
\alias{as.stratifiedvegdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reshapes community data from individual into stratified form
}
\description{
Function \code{stratifyvegdata} reshapes individual abundance values into species abundance values per size class or combination of size classes. Function \code{as.stratifiedvegdata} checks if the input list has appropriate properties and turns it into an object of class '\code{stratifiedvegdata}'.
}
\usage{
stratifyvegdata(x, sizes1, sizes2 = NULL, treeSel=NULL, spcodes = NULL, 
                plotColumn="plot", speciesColumn = "species", 
                abundanceColumn="abundance", size1Column = "size", size2Column = NULL, 
                cumulative=FALSE, counts=FALSE, mergeSpecies=FALSE, verbose=FALSE)
as.stratifiedvegdata(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data frame containing individual plant data. Individuals are in rows, while measurements are in columns.}
  \item{sizes1}{A numerical vector containing the breaks for primary size classes in ascending order.}
  \item{sizes2}{A numerical vector containing the breaks for secondary size classes in ascending order.}
  \item{treeSel}{A logical vector specifying which rows in \code{x} to be used. By default (\code{treeSel = NULL}) all rows are taken.}
  \item{spcodes}{A character vector indicating the codes of species to be used for stratification (species codes beyond those appearing in \code{x} are possible). If \code{spcodes = NULL} then all species in \code{x} are used.}
  \item{plotColumn}{The name of the column in \code{x} that contains plot identifiers.}
  \item{speciesColumn}{The name of the column in \code{x} that contains species names.}
  \item{abundanceColumn}{The name of the column in \code{x} that contains abundance values.}
  \item{size1Column}{The name of the column in \code{x} that contains values for primary size classes.}
  \item{size2Column}{The name of the column in \code{x} that contains values for secondary size classes.}
  \item{cumulative}{A flag to indicate that cumulative abundance profiles or surfaces are desired.}
  \item{counts}{A flag to indicate that the output should be individual counts instead of added abundance values.}
  \item{mergeSpecies}{A flag to indicate that species identity should be ignored. This leads to analyzing the structure of biomass disregarding species identity.}
  \item{verbose}{A logical flag to indicate extra output.}
  \item{X}{A list with as many elements as plot records. Each element should be of class 'matrix' or 'data.frame' with species in rows and strata in columns. Furthermore, the number of rows (species) and columns (strata) should be the same for all elements.}
}
\details{
For each individual (row) in \code{x}, \code{stratifyvegdata} assigns it to the size class (stratum) containing its size. The corresponding abundance value (e.g. crown cover) of the individual is added to the abundance of the corresponding species at the size class (stratum). If \code{sizes2} and \code{size2Column} are supplied, the function assigns each individual (row) in \code{x} to the combination of size classes (e.g. tree height and diameter). 
}
\value{
Both functions return an object of class '\code{stratifiedvegdata}', which is a list of matrices, one for each plot record. Each element (matrix) has as many rows as species and as many columns as size classes (i.e., as many as elements in vector \code{sizes1}). Columns are named starting with 'S' and continuing with the size class (stratum) number. If \code{mergeSpecies=TRUE} then all matrices have a single row (whose name is \code{"all"}). If \code{sizes2} and \code{size2Column} are supplied to \code{stratifyvegdata}, the function returns an object of class '\code{doublestratifiedvegdata}', which is a list of arrays, one for each plot record. Each element (array) has three dimensions corresponding to species, primary sizes (number of elements in in vector \code{sizes1}) and secondary sizes (number of elements in in vector \code{sizes2}). If \code{cumulative=TRUE} then the function returns cumulative abundances (see \code{\link{CAP}} and \code{\link{CAS}}).
}
\references{
De \enc{Cáceres}{Caceres}, M., Legendre, P. & He, F. (2013) Dissimilarity measurements and the size structure of ecological communities. Methods in Ecology and Evolution 4: 1167-1177.
}
\author{
Miquel De \enc{Cáceres}{Caceres}, CREAF.
}
\seealso{
\code{\link{reshape}}, \code{\link{CAP}}, \code{\link{CAS}}
}
\examples{
## Load tree data
data(treedata)

## Inspect tree data
head(treedata)

## Define stratum thresholds (4 strata)
heights = seq(0,4, by=0.5)
diameters = seq(0,2, by=0.5)

## Stratify tree data using heights as structural variable
X = stratifyvegdata(treedata, sizes1=heights, plotColumn="plotID",
                    speciesColumn="species", size1Column="height", counts=TRUE)

## Inspect the second plot record
X[[2]]

## Stratify tree data using heights as structural variable and cover as abundance
Y = stratifyvegdata(treedata, sizes1=heights, plotColumn="plotID",
                    speciesColumn="species", size1Column="height", 
                    abundanceColumn="cover")
Y[[2]]

## Stratify tree data using heights and diameters as structural variables
Z = stratifyvegdata(treedata, sizes1=heights, sizes2=diameters, plotColumn="plotID",
                    speciesColumn="species", size1Column="height", size2Column="diam",
                    counts=TRUE)
Z[[2]]



}
