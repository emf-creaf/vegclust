#' Reshapes as kmeans object
#'
#' This function casts an object of class \code{\link{vegclust}} into an object of class \code{\link{kmeans}}.
#' 
#' @param x An object of class \code{\link{vegclust}} to be casted, where \code{method="KM"} and \code{mode="raw"}.
#'
#' @returns An object of class \code{\link{kmeans}}
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' @seealso \code{\link{vegclust}}, \code{\link{kmeans}}
#' 
#' @export
#'
#' @examples
#' ## Loads data  
#' data(wetland)
#' 
#' ## This equals the chord transformation 
#' wetland.chord <- as.data.frame(sweep(as.matrix(wetland), 1, 
#'                                      sqrt(rowSums(as.matrix(wetland)^2)), "/"))
#' 
#' ## Create noise clustering with 3 clusters. Perform 10 starts from random seeds 
#' wetland.vc <- vegclust(wetland.chord, mobileCenters=3, 
#'                        method="KM", nstart=10)
#' 
#' ## Reshapes as kmeans object
#' wetland.km <- vegclust2kmeans(wetland.vc)
#' wetland.km
vegclust2kmeans<-function(x) {
  if(!inherits(x, "vegclust")) stop("x must be a vegclust object.")
  if(x$method!="KM" || x$mode!="raw") stop("Clustering model must be KM (kmeans) and mode must be raw.")
  cluster <- as.numeric(apply(x$memb==1,1,which))
  centers <- x$mobileCenters
  if(!is.null(x$fixedCenters)) centers<-cbind(x$fixedCenters,centers)
  rownames(centers) = 1:nrow(centers)
  withinss <- as.numeric(x$withinss)
  size <- as.numeric(x$size)
  iter <-1
  totss<-NULL
  ifault <- NULL
  tot.withinss<-sum(withinss)
  betweenss<-NULL
  cl<-list(cluster=cluster,centers=centers,totss=totss,
           withinss=withinss,tot.withinss=tot.withinss,betweenss=betweenss,
           size =size, iter=iter, ifault=ifault)
  class(cl)<-"kmeans"
  return(cl)
}