#' Turns into vegclust objects
#'
#' Attempts to turn its arguments into a \code{\link{vegclust}} object
#' 
#' @param x A site-by-species data matrix (raw mode), or a site-by-site distance matrix (distance mode).
#' @param y A vector indicating the cluster that each object in \code{x} belongs to. Alternatively, a fuzzy/hard site-by-group matrix of membership values.
#' @param method A clustering model from which \code{y} was obtained (normally "KM"). Current accepted models are:
#'   \itemize{
#'     \item{\code{"KM"}:}{ K-means or hard c-means (MacQueen 1967)}
#'     \item{\code{"KMdd"}:}{ Hard c-medoids (Krishnapuram et al. 1999)}
#'     \item{\code{"FCM"}:}{ Fuzzy c-means (Bezdek 1981)}
#'     \item{\code{"FCMdd"}:}{ Fuzzy c-medoids (Krishnapuram et al. 1999)}
#'     \item{\code{"NC"}:}{ Noise clustering (Dave and Krishnapuram 1997)}
#'     \item{\code{"NCdd"}:}{ Noise clustering with medoids}
#'     \item{\code{"HNC"}:}{ Hard noise clustering}
#'     \item{\code{"HNCdd"}:}{ Hard noise clustering with medoids}
#'     \item{\code{"PCM"}:}{ Possibilistic c-means (Krishnapuram and Keller 1993)}
#'     \item{\code{"PCMdd"}:}{ Possibilistic c-medoids}
#'   }
#' @param m The fuzziness exponent to be used, relevant for all fuzzy models (FCM, FCMdd, NC, NCdd, PCM and PCMdd).
#' @param dnoise The distance to the noise cluster, relevant for noise clustering models (NC, HNC, NCdd and HNCdd). 
#' @param eta A vector of reference distances, relevant for possibilistic models (PCM and PCMdd).
#'
#' @details
#' This function is used to generate \code{\link{vegclust}} objects which can then be used in \code{\link{vegclass}} to classify new data. If the input classification is hard (i.e. yes/no membership), cluster centers are calculated as multivariate means, and the method for assigning new data is assumed to be k-means (\code{"KM"}), i.e. plots will be assigned to the nearest cluster center. If community data is given as site-by-species data matrix the cluster centroids are added as \code{mobileCenters} in the \code{\link{vegclust}} object. Centroids will not be computed if community data is given as a site-by-site dissimilarity matrix. Moreover, current implementation does not allow \code{y} to be a membership matrix when \code{x} is a distance matrix.
#' 
#' @returns An object of class \code{\link{vegclust}}.
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF.
#' @seealso \code{\link{vegclust}}, \code{\link{vegclass}}, \code{\link[vegan]{decostand}}
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
#' ## Splits wetland data into two matrices of 30x27 and 11x22
#' wetland.30 <- wetland.chord[1:30,]
#' wetland.30 <- wetland.30[,colSums(wetland.30)>0]
#' dim(wetland.30)
#' wetland.11 <- wetland.chord[31:41,]
#' wetland.11 <- wetland.11[,colSums(wetland.11)>0] 
#' dim(wetland.11)
#' 
#' ## Performs a K-means clustering of the data set with 30 sites
#' wetland.km <- kmeans(wetland.30, centers=3, nstart=10)
#' 
#' ## Transforms the 'external' classification of 30 sites into a 'vegclust' object
#' wetland.30.vc <- as.vegclust(wetland.30, wetland.km$cluster)
#' 
#' ## Assigns the second set of sites according to the (k-means) membership rule 
#' ## That is, sites are assigned to the cluster whose cluster centroids is nearest.
#' wetland.11.km <- vegclass(wetland.30.vc, wetland.11)
#' 
#' ## A similar 'vegclust' object is obtained when using the distance mode...
#' wetland.d.vc <- as.vegclust(dist(wetland.30), wetland.km$cluster)
#' 
#' ## which can be also used to produce the assignment of the second set of objects
#' wetland.d.11 <- as.data.frame(as.matrix(dist(wetland.chord)))[31:41,1:30]
#' wetland.d.11.km <- vegclass(wetland.d.vc,wetland.d.11)
as.vegclust <-
function(x,y, method="KM", m=1.0, dnoise=NULL, eta=NULL) {
  METHODS <- c("KM", "FCM", "PCM","NC","HNC" ,"KMdd","NCdd", "HNCdd", "FCMdd", "PCMdd")
  method <- match.arg(method, METHODS)
  
   dist2onecluster<-function(x,object) {
		#x is an (euclidean) distance matrix
		x = as.matrix(x)
		n = nrow(x)
		if (length(object)!=n) stop("Length of object vector must be equal to the number of sites in x")
		vargeom = (rep(1,n) %*% (x^2) %*% rep(1,n))/(2*(n^2))
		return(sqrt((sum(object^2)/n)-vargeom))
	}
   dist2clusters<-function(x,cluster,object) {
	   n = nrow(as.matrix(x))
		if (length(cluster)!=n) 
            stop("Length of cluster vector must be equal to the number of sites in x")
		cluster = as.factor(cluster)
		k = length(levels(cluster))
		d = vector("numeric",k)
		for(i in 1:k) {
			sel = (cluster==levels(cluster)[i])
		  sel[is.na(sel)]=FALSE
			d[i] = dist2onecluster(as.dist(as.matrix(x)[sel,sel]),object[sel])
		}	
		names(d) = levels(cluster)
		return(d)
   }	
	
   if(inherits(x,"dist")) {
   	mode="distance"
   	x = as.matrix(x)
   	sitenames = rownames(x)
   } else {
   	mode="raw"
   	sitenames = rownames(x)
   	varnames = names(x)
   	x = as.matrix(x)
   }
   if(is.vector(y)) {
   	cluster = y
   	cln =levels(as.factor(cluster))
    u = as.memb(cluster)
   	rownames(u) = sitenames
    colnames(u) = cln
   	k = length(cln)
   	if(method=="NC"||method=="NCdd"|| method=="HNC"||method=="HNCdd") {
       u = cbind(u, rep(0,nrow(u)))
       colnames(u)[k+1] = "N"
   	}
    u = as.data.frame(u)
   } else if(is.matrix(y) || is.data.frame(y)) {
   	u = as.data.frame(y)
   	cln = names(u)
   	k = length(cln)
   	if(method=="NC"||method=="NCdd"|| method=="HNC"||method=="HNCdd") {
   	  k = k-1
   	  cln = cln[1:k]
   	}
   }
   n = nrow(x)
   
   dist2cent = matrix(0,nrow=n,ncol=k) 
      
   if(mode=="distance") {
   	for(j in 1:n) {
   		dist2cent[j,] = dist2clusters(x,cluster,x[j,])
   	}
   	centers=NULL
   } else if (mode=="raw"){
   	cm = clustcentroid(x,u[,1:k], m=m)
   	colnames(cm) = varnames
   	rownames(cm) = cln
   	centers = as.data.frame(cm)
   	for(i in 1:k) {
   		dist2cent[,i] = sqrt(rowSums(sweep(x,2,cm[i,],"-")^2))
   	}   
   }   
   dist2cent = as.data.frame(dist2cent)  
   names(dist2cent) = cln
   rownames(dist2cent) = sitenames
	
   size = colSums(u[,1:k])
   withinss = colSums((dist2cent^2)*u[,1:k])
   functional = sum(withinss)
	
   res = list(mode = mode, method=method, m = m, dnoise = dnoise, eta = eta, memb=u,mobileCenters=centers, fixedCenters=NULL, dist2clusters=dist2cent, withinss = withinss, size=size, functional=functional)
   class(res)<-"vegclust"
	return(res)
}

