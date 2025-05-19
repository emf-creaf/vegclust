#' Cross-table of two fuzzy classifications
#' 
#' Calculates a cross-tabulated matrix relating two fuzzy membership matrices
#'
#' @param x A site-by-group fuzzy membership matrix. Alternatively, an object of class 'vegclust' or 'vegclass'.
#' @param y A site-by-group fuzzy membership matrix. Alternatively, an object of class 'vegclust' or 'vegclass'.
#' @param relativize If \code{TRUE} expresses the cross-tabulated values as proportions of cluster size in \code{x}.
#'
#' @returns
#' A cross-tabulated matrix comparing the two classifications. In general, each cell's value is the (fuzzy) number of objects that in \code{x} are assigned to the cluster corresponding to the row and in \code{y} are assigned to the cluster corresponding to the column. If \code{relativize=TRUE} then the values of each row are divided by the (fuzzy) size of the corresponding cluster in \code{x}. 
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF.
#' 
#' @seealso \code{\link{defuzzify}}, \code{\link{vegclust}}, \code{\link[vegan]{decostand}}
#' @export
#'
#' @examples
#' ## Loads data  
#' data(wetland)
#'   
#' ## This equals the chord transformation 
#' wetland.chord <- as.data.frame(sweep(as.matrix(wetland), 1, 
#'                                sqrt(rowSums(as.matrix(wetland)^2)), "/"))
#' 
#' ## Create clustering with 3 clusters. Perform 10 starts from random seeds 
#' ## and keep the best solution. Try both FCM and NC methods:
#' wetland.fcm <- vegclust(wetland.chord, mobileCenters=3, m = 1.2, method="FCM", nstart=10)
#' wetland.nc <- vegclust(wetland.chord, mobileCenters=3, m = 1.2, dnoise=0.75, method="NC", 
#'                       nstart=10)
#' 
#' ## Compare the results
#' crossmemb(wetland.fcm, wetland.nc, relativize=FALSE)
crossmemb<-function(x,y,relativize=TRUE) {
	if(inherits(x, "vegclust") || inherits(x, "vegclass")) x = x$memb
	if(inherits(y, "vegclust") || inherits(y, "vegclass")) y = y$memb
	c=t(x)%*%as.matrix(y)
	if(relativize) c = sweep(c,1,colSums(x),"/")
    return(c)
}