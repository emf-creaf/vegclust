#' Distance between pairs of cluster centroids
#'
#' Calculates the distance between pairs of cluster centroids, given a distance matrix and a cluster vector.
#' 
#' @param x A site-by-site data matrix or an object of class \code{\link{dist}} containing the distance values between pairs of sites (plot records).
#' @param cluster A vector indicating the hard membership of each object in \code{x} to a set of groups. Can contain \code{NA} values.
#'
#' @returns
#' An object of class \code{\link{dist}} containing the distances between cluster centers.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' 
#' @export
#'
interclustdist <-
function(x,cluster) {
	  D = as.dist(x)
	  memb = as.memb(cluster)
	  vars<-clustvar(D,cluster)

	  disp1<-function(f1,f2, D, var1,var2) {
        if (is.na(sum(f1)) || is.na(sum(f2))) cd <- NA
        else if (sum(f1) < 1e-16 || sum(f2) < 1e-16) cd <- 0
        else {
        		cd <- sqrt(((f1 %*% (as.matrix(D)^2) %*% f2)/(sum(f1)*sum(f2))) - var1-var2)
    		 }
    		 return(cd)
    	}
    	k <- ncol(memb)
    	interD = matrix(0,nrow=k,ncol=k)
    	for(i in 1:(k-1)) {
	    	for(j in (i+1):k) {
  	  			interD[i,j] = interD[j,i] = disp1(memb[,i], memb[,j],D, vars[i],vars[j])
    		}
    	}
    	return(as.dist(interD))
}