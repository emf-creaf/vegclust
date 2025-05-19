#' Fuzzy Clustering of Vegetation Data
#'
#' A set of functions to: (1) perform fuzzy clustering of vegetation data; (2) to assess ecological community similarity on the basis of structure and composition.
#'
#' @name vegclust-package
#' @aliases vegclust-package
#' @docType package
#' @author \strong{Maintainer}: Miquel De Cáceres
#' \email{miquelcaceres@@gmail.com}
#' [ORCID](https://orcid.org/0000-0001-7132-2080)
#'
#' @seealso Useful links: \itemize{ \item{
#' \url{https://emf-creaf.github.io/vegclust/index.html}} }
#'
#' @references De Caceres et al, 2010 (\doi{10.1111/j.1654-1103.2010.01211.x}), De Caceres et al, 2013 (\doi{10.1111/2041-210X.12116}).
#' @keywords internal
#' @examples
#' ## Loads data  
#' data(wetland)
#' 
#' ## This equals the chord transformation 
#' wetland.chord <- as.data.frame(sweep(as.matrix(wetland), 1, 
#'                                      sqrt(rowSums(as.matrix(wetland)^2)), "/"))
#' 
#' ## Create noise clustering with 3 clusters. Perform 10 starts from random seeds 
#' ## and keep the best solution
#' wetland.nc <- vegclust(wetland.chord, mobileCenters=3, m = 1.2, dnoise=0.75, 
#'                        method="NC", nstart=10)
#' 
"_PACKAGE"

## usethis namespace: start
#' @importFrom grDevices colorRampPalette rainbow
#' @importFrom graphics abline axis legend lines matlines persp arrows text
#' @importFrom vegan decostand
#' @importFrom stats as.dist cutree dist quantile model.matrix var cmdscale
## usethis namespace: end
NULL