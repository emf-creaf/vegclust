#' Noise clustering with increasing number of clusters
#' 
#' Performs several runs of function 'vegclust' on a community data matrix using an increasing number of clusters until some conditions are met.
#'
#' @param x Community data table. A site (rows) by species (columns) matrix or data frame.
#' @param method A clustering model. Current accepted models are of the noise clustering family: 
#' \itemize{
#'   \item{\code{"NC"}: Noise clustering (Dave and Krishnapuram 1997)}
#'   \item{\code{"NCdd"}: Noise clustering with medoids}
#'   \item{\code{"HNC"}: Hard noise clustering}
#'   \item{\code{"HNCdd"}: Hard noise clustering with medoids}
#' }  
#' @param ini.fixed.centers The coordinates of initial fixed cluster centers. These will be used as \code{fixedCenters} in all calls to \code{\link{vegclust}}. If \code{method="NCdd"} or \code{method="HNCdd"} then \code{ini.fixed.centers} can be specified as a vector of indices for medoids.
#' @param min.size The minimum size (cardinality) of clusters. If any of the current k clusters does not have enough members the algorithm will stop and return the solution with k-1 clusters.
#' @param max.var The maximum variance allowed for clusters (see function \code{\link{clustvar}}). If specified, the algorithm will stop when any of the clusters is at the same time small and has large variance. If \code{max.var = NULL} then this criterion is not used.
#' @param alpha Criterion to choose cluster seeds from the noise class. Specifically, an object is considered as cluster seed if the membership to the noise class is larger than \code{alpha}.
#' @param nstart A number indicating how many random trials should be performed for number of groups. Each random trial uses the k-1 cluster centers plus the coordinates of the current cluster seed as initial solution for \code{\link{vegclust}}. Thus, if there are less cluster seed candidates than \code{nstart}, then not all runs are conducted.
#' @param fix.previous Flag used to indicate that the cluster centers found when determining k-1 clusters are fixed when determining k clusters.
#' @param dnoise The distance to the noise cluster. 
#' @param m The fuzziness exponent.
#' @param ... Additional parameters for function \code{\link{vegclust}}.
#'
#' @details
#' Function \code{hier.vegclust} takes starting cluster configurations from cuts of a dendrogram given by object \code{hclust}. Function \code{random.vegclust} chooses random objects as cluster centroids and for each number of clusters performs \code{nstart} trials.
#' 
#' @returns
#' Returns an object of class \code{\link{vegclust}}; or \code{NULL} if the initial cluster does not contain enough members.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' 
#' @references
#' \enc{Davé}{Dave}, R. N. and R. Krishnapuram (1997) Robust clustering methods: a unified view. IEEE Transactions on Fuzzy Systems 5, 270-293.
#' 
#' @seealso \code{\link{vegclust}},\code{\link{hier.vegclust}} 
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
#' ## Call incremental noise clustering 
#' wetland.nc <- incr.vegclust(wetland.chord, method="NC", m = 1.2, dnoise=0.75, 
#'                             min.size=5)
#' 
#' ## Inspect cluster sizes
#' print(wetland.nc$size)
incr.vegclust<-function(x, method="NC", ini.fixed.centers = NULL, 
                        min.size = 10, max.var = NULL, alpha = 0.5, 
                        nstart=100, fix.previous = TRUE, dnoise=0.75, m=1.0,...) {
  METHODS <- c("NC", "HNC", "NCdd", "HNCdd")
  method <- match.arg(method, METHODS)  
  k=1
  if(is.null(ini.fixed.centers)) {
    cat(paste("Vegclust with one new group..."))
    vc<-vegclust(x,mobileCenters=k, method=method, nstart=nstart,dnoise=dnoise,m=m,...)
  } else { #Assign plots and look for unassigned ones
    y = as.matrix(x)
    if(is.vector(ini.fixed.centers) && (method=="HNCdd" || method=="NCdd")) ini.fixed.centers = x[ini.fixed.centers,]
    else ini.fixed.centers = as.matrix(ini.fixed.centers)
    dist2cent = matrix(0,nrow=nrow(y),ncol=nrow(ini.fixed.centers))
    for(i in 1:nrow(ini.fixed.centers)) {
      dist2cent[,i] = sqrt(rowSums(sweep(y,2,ini.fixed.centers[i,],"-")^2))
    }
    if(method=="NC" || method=="NCdd") {
       d2cm2<-cbind(dist2cent,dnoise)^2
       a<-sweep(d2cm2,1,d2cm2[,ncol(d2cm2)],"/")
       seeds<-which((1/rowSums(a^(-1/(m-1))))>alpha)
    } else if(method=="HNC" || method=="HNCdd") {
      d2cm<-cbind(dist2cent,dnoise)
      minC<-apply(d2cm,1,which.min)
      seeds<-which(minC==ncol(d2cm))
    }
    cat(paste("Number of initial groups: ", nrow(ini.fixed.centers),"\n", sep=""))
    cat(paste("Number of initial seeds: ", length(seeds),"\n", sep=""))
    if(nstart<length(seeds)) seeds = sample(seeds,nstart)
    cat(paste("Vegclust with 1 new group"))
    vcbest = NULL
    for(i in 1:length(seeds)) {
      cat(".")
      vc<-vegclust(x,mobileCenters=x[seeds[i],],
                   fixedCenters=ini.fixed.centers, 
                   method=method,dnoise=dnoise, m=m,...)
      if(is.null(vcbest)) vcbest = vc
      else if(vc$functional<vcbest$functional) vcbest = vc
    }
    vc = vcbest
  }
  #Number of seed objects 
  noise<-vc$memb[,ncol(vc$memb)]>alpha
  cat(paste("Number of remaining cluster seeds: ", sum(noise),"\n",sep=""))
  #Stop before continuing if there are not enough seeds object or the first cluster is too small
  if(sum(colSums(vc$memb[,-ncol(vc$memb), drop=FALSE]>alpha)<min.size)>0) {
    cat("The initial cluster was too small.\n")
    return(NULL)
  }
  if(sum(noise)==0) {
    cat("Not enough objects to act as cluster seeds. Stopping.\n")
    return(vc)
  } 
  #Continue
  cont = TRUE
	while(cont) {
	  k = k+1
	  vcold = vc
	  vcbest = NULL
    if(nstart<sum(noise)) seeds = sample(which(noise),nstart)
    else seeds = which(noise)
	  cat(paste("Vegclust with", k,"new groups"))
    fixed = vcold$mobileCenters[1,]
	  if(!is.null(vcold$fixedCenters)) fixed = rbind(fixed,vcold$fixedCenters)
	  for(i in 1:length(seeds)) {
	    cat(".")
	    if(fix.previous) vc<-vegclust(x,mobileCenters=x[seeds[i],],
                                    fixedCenters=fixed, method=method, dnoise=dnoise, m=m,...)
      else vc<-vegclust(x,mobileCenters=rbind(fixed,x[seeds[i],]), 
                        fixedCenters = ini.fixed.centers, method=method, dnoise=dnoise, m=m,...)
      if(is.null(vcbest)) vcbest = vc
      else if(vc$functional<vcbest$functional) vcbest = vc
    }
    vc = vcbest
	  noise<-vc$memb[,ncol(vc$memb)]>alpha
	  cat(paste("Number of remaining cluster seeds: ", sum(noise),"\n",sep=""))    
    nboth=0
    nsmall=0
    if(!is.null(max.var)) {
       nboth = sum(clustvar(vc)>max.var & colSums(vc$memb[,-ncol(vc$memb), drop=FALSE]>alpha)<min.size)
       cont = (nboth==0 && sum(noise)>0)
    } else {
      nsmall = sum(colSums(vc$memb[,-ncol(vc$memb), drop=FALSE]>alpha)<min.size)
      cont = (nsmall==0 && sum(noise)>0)
    }
	}
	if(is.null(max.var) && nsmall>0) {
	  cat("Some of the current clusters are too small. ")
	  cat(paste("Returning vegclust with", k-1,"new group(s).\n"))
	  return(vcold)
	} else if(!is.null(max.var) && nboth>0) {
	    cat("Some of the current clusters are too small and have large variance. ")
	    cat(paste("Returning vegclust with", k-1,"new group(s).\n"))
	    return(vcold)
	} else if(sum(noise)==0) {
	  cat("Not enough objects to act as cluster seeds. Stopping.\n")
	  return(vc)
	}
  return(NULL)
}