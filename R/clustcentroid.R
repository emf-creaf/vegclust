#' Cluster centers of a classification
#' 
#' Function \code{clustcentroid} calculates the centroid (multivariate average) coordinates of a classification. Function \code{clustmedoid} determines the medoid (object whose average dissimilarity to all the other objects is minimal) for each cluster in the classification.
#'
#' @param x Community data, a site-by-species data frame. In function \code{clustmedoid}, \code{x} can alternatively be an object of class \code{\link{dist}} (otherwise, the dissimilarity measure is assumed to be the Euclidean distance).
#' @param y It can be (a) A vector indicating the cluster that each object in \code{x} belongs to; (b) a fuzzy/hard site-by-group matrix of membership values; (c) an object of class \code{\link{vegclust}} or \code{\link{vegclass}}
#' @param m Fuzziness exponent, only effective when \code{y} is a fuzzy membership matrix.
#'
#' @returns
#' Function \code{clustcentroid} returns a group-by-species matrix containing species average abundance values (i.e. the coordinates of each cluster centroid). Function \code{clustmedoid} returns a vector of indices (medoids).
#' 
#' @details
#' In order to assign new plot record data into a predefined set of classes, one should use functions \code{\link{as.vegclust}} and \code{\link{vegclass}} instead.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF.
#' 
#' @seealso \code{\link{as.vegclust}}, \code{\link{vegclass}}, \code{\link{vegclust}}, \code{\link{kmeans}}
#' 
#' @export
#'
#' @name clustcentroid
#' @examples
#' ## Loads stats
#' library(stats)
#' 
#' ## Loads data
#' data(wetland)
#' 
#' ## This equals the chord transformation 
#' ## (see also \code{\link{decostand}} in package 'vegan')
#' wetland.chord <- as.data.frame(sweep(as.matrix(wetland), 1,
#'                                      sqrt(rowSums(as.matrix(wetland)^2)), "/"))
#' 
#' ## Performs a K-means clustering
#' wetland.km <- kmeans(wetland.chord, centers=3, nstart=10)
#' 
#' ## Gets the coordinates corresponding to the centroids of KM clusters
#' clustcentroid(wetland.chord, y=wetland.km$cluster)
#' 
#' ## Gets the object indices corresponding to the medoids of KM clusters
#' clustmedoid(wetland.chord, y=wetland.km$cluster)
clustcentroid <-
function(x,y, m=1) {
  if(inherits(y,"vegclust") || inherits(y, "vegclass")) {
    u = as.matrix(y$memb)^y$m
    if(y$method=="NC"||y$method=="HNC"||y$method=="NCdd"||y$method=="HNCdd") {
      u = u[,-ncol(u)]
    }
  } else if(is.null(dim(y))) {
	  u = as.memb(y)
	  colnames(u) = levels(as.factor(y))
  } else {
		u = as.matrix(y)^m
	}
  s = t(as.matrix(x))%*%(u)
  centers = sweep(t(s),1,colSums(u),"/")
  colnames(centers) = names(x)
  rownames(centers) = colnames(u)   
	return(centers)
}

#' @rdname clustcentroid
#' @export 
clustmedoid <-
  function(x, y, m=1) {
    if(inherits(y,"vegclust") || inherits(y, "vegclass")) {
      u = as.matrix(y$memb)^y$m
      if(y$method=="NC"||y$method=="HNC"||y$method=="NCdd"||y$method=="HNCdd") {
        u = u[,-ncol(u)]
      }
    } else if(is.null(dim(y))) {
      u = as.memb(y)
      colnames(u) = levels(as.factor(y))
    } else {
      u = as.matrix(y)^m
    }
    c = ncol(u)
    med = numeric(c)
    if(!inherits(x,"dist")) {
      d = as.matrix(dist(x))    
    }
    for(k in 1:c) {
      med[k] = which.min((u[,k]^m)%*%d)
    }   
    names(med) = row.names(as.matrix(d))[med]
    return(med)
  }
