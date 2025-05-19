#' Defuzzifies a fuzzy partition
#' 
#' Transforms a fuzzy classification into a crisp (hard) classification.
#'
#' @param object A site-by-group fuzzy membership matrix. Alternatively, an object of class 'vegclust' or 'vegclass'.
#' @param method Either \code{"max"} to choose for the maximum membership value across clusters, or \code{"cut"} for an alpha-cut.
#' @param alpha Threshold for the alpha-cut, bounded between 0 and 1.
#' @param na.rm If \code{TRUE} removes the objects that do not belong to any cluster when using \code{method="cut"}.
#'
#' @details
#' Alpha-cut means that memberships lower than alpha are transformed into 0 while memberships higher than alpha are transformed into 1. This means that if alpha values are low (i.e. close to 0), an object may belong to more than one group after defuzzification. These will generate a concatenation of cluster names in the output \code{cluster} vector and a row with sum more than one in the \code{memb} matrix). Similarly, if alpha is high (i.e. close to 1) there are objects that may be left unclassified. These will get \code{NA} in the \code{cluster} vector and zero row in the \code{memb} matrix.
#' 
#' @returns
#' A list with the following items:
#' \itemize{
#'   \item{\code{memb}: A data frame with the hard membership partition.}
#'   \item{\code{cluster}: A vector (factor) with the name of the cluster for each object.}
#' }
#' 
#' @references
#' \enc{Davé}{Dave}, R. N. and R. Krishnapuram (1997) Robust clustering methods: a unified view. IEEE Transactions on Fuzzy Systems 5, 270-293.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF.
#' 
#' @seealso \code{\link{vegclust}}, \code{\link[vegan]{decostand}}
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
#' ## and keep the best solution
#' wetland.nc <- vegclust(wetland.chord, mobileCenters=3, m = 1.2, dnoise=0.75, 
#'                        method="NC", nstart=10)
#' 
#' ## Defuzzification using an alpha-cut (alpha=0.5)
#' wetland.nc.df <- defuzzify(wetland.nc$memb, method="cut")
#' 
#' ## Cluster vector, with 'N' for objects that are unclassified, 
#' ## and 'NA' for objects that are intermediate
#' print(wetland.nc.df$cluster)
#' 
#' ## Hard membership matrix (site 22 does not get any cluster assigned)
#' print(wetland.nc.df$memb)
defuzzify<-function(object, method="max", alpha=0.5,na.rm=FALSE) {
    METHODS <- c("max", "cut")
    method <- match.arg(method, METHODS)	
	if(inherits(object,"vegclust")|| inherits(object,"vegclass")) memb = object$memb
	else memb = object
	
	memb = as.data.frame(memb)
	clnames = names(memb)
	cluster = vector("character",nrow(memb))
	u = matrix(0,nrow(memb), ncol(memb))
	cluster = rep(NA,nrow(memb))
	
	a<-function(v) {
		return(paste(clnames[which(v==1)],collapse="+"))
	}
  if(method=="max") {
  	  vmax<-apply(memb,1,max)
   	u=ifelse(memb==vmax,1,0)
  } else if (method=="cut") {
   	u[memb>alpha]= 1
  }
  cluster = as.character(apply(u,1,a))
  cluster[cluster==""]= NA

	u = as.data.frame(u)
	row.names(u) = row.names(memb)
	names(u) = names(memb)
	if(na.rm) {
		sel = !is.na(cluster)
		cluster =cluster[sel]
		u = subset(u,sel)
	}
  names(cluster) = row.names(u)
	return(list(memb=u,cluster=cluster))
}