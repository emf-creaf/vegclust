\encoding{UTF-8}
\name{vegclass}
\alias{vegclass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Classifies vegetation communities 
}
\description{
Classifies vegetation communities into a previous fuzzy or hard classification.
}
\usage{
vegclass(y, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
	An object of class \code{\link{vegclust}} that represents a previous knowledge.
   }
  \item{x}{
Community data to be classified, in form of a site by species matrix (if the vegclust object is in \code{raw} mode) or a data frame containing the distances between the new sites in rows and the old sites in columns (if the vegclust object is in \code{distance} mode).
   }
}
\details{
This function uses the classification model specified in \code{y} to classify the communities (rows) in \code{x}. When vegclust is in \code{raw} mode, the function calls first to \code{\link{conformveg}} in order to cope with different sets of species. See the help of \code{\link{as.vegclust}} to see an example of \code{vegclass} with distance matrices.
}
\value{
Returns an object of type \code{vegclass} with the following items:
  \item{method}{The clustering model used in \code{y}}
  \item{m}{The fuzziness exponent in \code{y}}
  \item{dnoise}{The distance to the noise cluster used for noise clustering (models NC, NCdd, HNC, HNCdd). This is set to \code{NULL} for other models.}
  \item{eta}{ The reference distance vector used for possibilistic clustering (models PCM and PCMdd). This is set to \code{NULL} for other models.}
  \item{memb}{The fuzzy membership matrix.}
  \item{dist2clusters}{The matrix of object distances to cluster centers.}
}
\references{
\enc{Davé}{Dave}, R. N. and R. Krishnapuram (1997) Robust clustering methods: a unified view. IEEE Transactions on Fuzzy Systems 5, 270-293.

Bezdek, J. C. (1981) Pattern recognition with fuzzy objective functions. Plenum Press, New York.

Krishnapuram, R. and J. M. Keller. (1993) A possibilistic approach to clustering. IEEE transactions on fuzzy systems 1, 98-110.

De \enc{Cáceres}{Caceres}, M., Font, X, Oliva, F. (2010) The management of numerical vegetation classifications with fuzzy clustering methods [Related software]. Journal of Vegetation Science 21 (6): 1138-1151.
}
\author{
Miquel De \enc{Cáceres}{Caceres}, Forest Science Center of Catalonia.
}
\seealso{
\code{\link{vegclust}}, \code{\link{as.vegclust}}, \code{\link{kmeans}}, \code{\link[e1071]{cmeans}}, \code{\link{conformveg}}
}
\examples{
## Loads data (38 columns and 33 species)
data(wetland)
dim(wetland)
  
## This equals the chord transformation 
## (see also  \code{\link{decostand}} in package vegan)
wetland.chord = as.data.frame(sweep(as.matrix(wetland), 1, 
                              sqrt(rowSums(as.matrix(wetland)^2)), "/"))

## Splits wetland data into two matrices of 30x27 and 11x22
wetland.30 = wetland.chord[1:30,]
wetland.30 = wetland.30[,colSums(wetland.30)>0]
dim(wetland.30)
wetland.11 = wetland.chord[31:41,]
wetland.11 = wetland.11[,colSums(wetland.11)>0] 
dim(wetland.11)

## Create noise clustering with 3 clusters from the data set with 30 sites. 
wetland.30.nc = vegclust(wetland.30, mobileCenters=3, m = 1.2, dnoise=0.75,
                         method="NC", nstart=10)

## Cardinality of fuzzy clusters (i.e., the number of objects belonging to)
wetland.30.nc$size

## Classifies the second set of sites according to the clustering of the first set
wetland.11.nc = vegclass(wetland.30.nc, wetland.11)

## Fuzzy membership matrix
wetland.11.nc$memb

## Obtains hard membership vector, with 'N' for objects that are unclassified
defuzzify(wetland.11.nc$memb)$cluster

}
