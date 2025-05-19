#' Fuzzy evaluation statistics
#' 
#' Computes several evaluation statistics on the fuzzy clustering results on objects of class \code{\link{vegclust}}. 
#'
#' @param y An object of class \code{\link{vegclust}} or a membership matrix.
#'
#' @details
#' These statistics were conceived to be computed on fuzzy partitions, such as the ones coming from Fuzzy C-means (Bezdek 1981). Maximum values of PCN or minimum values of PEN can be used as criteria to choose the number of clusters.
#' 
#' @returns
#' Returns an vector of four values: partition coefficient (PC), normalized partition coefficient (PCN), partition entropy (PE) and normalized partition entropy (PEN).
#' 
#' @references
#' Bezdek, J. C. (1981) Pattern recognition with fuzzy objective functions. Plenum Press, New York.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF.
#' 
#' @seealso \code{\link{vegclust}}
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
#' ## Create noise clustering with 2, 3 and 4 clusters. Perform 10 starts from random seeds 
#' ## and keep the best solutions
#' wetland.fcm2 <- vegclust(wetland.chord, mobileCenters=2, m = 1.2, method="FCM", nstart=10)
#' wetland.fcm3 <- vegclust(wetland.chord, mobileCenters=3, m = 1.2, method="FCM", nstart=10)
#' wetland.fcm4 <- vegclust(wetland.chord, mobileCenters=4, m = 1.2, method="FCM", nstart=10)
#' 
#' ## Compute statistics. Both PCN and PEN indicate that three groups are more advisable 
#' ## than 2 or 4.
#' print(vegclustIndex(wetland.fcm2))
#' print(vegclustIndex(wetland.fcm3))
#' print(vegclustIndex(wetland.fcm4))
vegclustIndex <-function (y) {
	
    partition.coefficient <- function(clres) {
    	   if(inherits(clres,"vegclust")) clres = clres$memb
        xrows <- dim(clres)[1]
        partitioncoefficient <- sum(apply(clres^2, 1, sum))/xrows
        return(partitioncoefficient)
    }
    partition.entropy <- function(clres) {
    	   if(inherits(clres,"vegclust")) clres = clres$memb
        xrows <- dim(clres)[1]
        ncenters <- dim(clres)[2]
        partitionentropy <- 0
        for (i in 1:xrows) {
            for (k in 1:ncenters) {
                if (clres[i, k] != 0) 
                  partitionentropy <- partitionentropy + (clres[i,k] * log(clres[i, k]))
            }
        }
        partitionentropy <- partitionentropy/((-1) * xrows)
        return(partitionentropy)
    }
    if(inherits(y,"vegclust")) y = y$memb

    xrows <- dim(y)[1]
    ncenters <- dim(y)[2]
    findex = numeric(4)
    findex[1] = partition.coefficient(y)
    findex[2] = (findex[1]*ncenters-1)/(ncenters-1)
    findex[3] = partition.entropy(y)
    findex[4] = findex[3]/(1-(ncenters/xrows))
    names(findex)<-c("PC", "PCN","PE","PEN")
    return(findex)
}


