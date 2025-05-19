#' Clustering with several number of clusters
#' 
#' Performs several runs of function 'vegclust' (or 'vegclustdist') on a community data matrix (or distance matrix) using different number of clusters
#'
#' @param x For \code{hier.vegclust} and \code{random.vegclust}, a site (rows) by species (columns) matrix or data frame. For \code{hier.vegclustdist} and \code{random.vegclustdist}, a square distance matrix.
#' @param hclust A hierarchical clustering represented in an object of type \code{\link{hclust}}.
#' @param cmin Number of minimum mobile clusters.
#' @param cmax Number of maximum mobile clusters.
#' @param min.size If \code{min.size != NULL}, it specifies the minimum size of clusters. If some clusters are smaller, the algorithm will return the solutions corresponding to lower numbers of clusters.
#' @param verbose Flag used to print which number of clusters is currently running.
#' @param ... Additional parameters for function \code{\link{vegclust}} or \code{\link{vegclustdist}}.
#'
#' @details
#' Function \code{hier.vegclust} takes starting cluster configurations from cuts of a dendrogram given by object \code{hclust}. Function \code{random.vegclust} chooses random objects as cluster centroids and for each number of clusters performs \code{nstart} trials. Functions \code{hier.vegclustdist} and \code{random.vegclustdist} are analogous to \code{hier.vegclust} and \code{random.vegclust} but accept distance matrices as input.
#' 
#' @returns
#' Returns an object of type 'mvegclust' (multiple vegclust), which contains a list vector with objects of type \code{\link{vegclust}}.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' 
#' @seealso \code{\link{vegclust}}, \code{\link{vegclustdist}}, \code{\link{vegclass}}, \code{\link{defuzzify}}, \code{\link{hclust}}
#' 
#' @export
#'
#' @examples
#' ## Loads data  
#' data(wetland)
#' 
#' ## This equals the chord transformation 
#' wetland.chord <- as.data.frame(sweep(as.matrix(wetland), 1, 
#'                                     sqrt(rowSums(as.matrix(wetland)^2)), "/"))
#' 
#' ## Create noise clustering from hierarchical clustering at different number of clusters
#' wetland.hc <- hclust(dist(wetland.chord),method="ward") 
#' wetland.nc1 <- hier.vegclust(wetland.chord, wetland.hc, cmin=2, cmax=5, 
#'                              m = 1.2, dnoise=0.75, method="NC")
#' 
#' ## Create noise clustering from random seeds at different levels
#' wetland.nc2 <- random.vegclust(wetland.chord, cmin=2, cmax=5, nstart=10,
#'                                m = 1.2, dnoise=0.75, method="NC")
hier.vegclust<-function(x, hclust, cmin=2,cmax=20, min.size = NULL,  verbose=TRUE,...) {
  vc = vector("list",length=cmax-cmin +1)
  xs = x[row.names(x)%in% hclust$labels,] #Select those objects used for hierarchical clustering
  i=1
  for(c in cmin:cmax) {
    if(verbose) cat(paste("PROCESSING",c,"MOBILE CLUSTERS\n"))
    if(c==1) {
      cent=matrix(0,nrow=1, ncol=ncol(xs))
      cent[1,]=colMeans(xs)
    }
    else cent = as.vegclust(xs,cutree(hclust,k=c))$mobileCenters
    vc[[i]] = vegclust(x,mobileCenters=cent,...)
    if(!is.null(min.size)) {
      nsmall = sum(vc[[i]]$size< min.size)
      if(nsmall>0) {
        cat("At least one cluster was too small. Stopping.")
        if(i==1) vc<-vector("list",length=0)
        else vc<-vc[1:(i-1)]
        break
      }
    }
    i = i+1
  }
  class(vc)<-"mvegclust"
  return(vc)
}


#' @rdname hier.vegclust
#' @param nstart A number indicating how many random trials should be performed for each number of groups
#' 
#' @export
random.vegclust<-function(x, cmin=2, cmax=20, nstart=10,  min.size = NULL, verbose=TRUE,...) {
  nvc = cmax-cmin +1
  vc = vector("list",nvc)
  i=1
  for(c in cmin:cmax) {
    if(verbose) cat(paste("PROCESSING",c,"MOBILE CLUSTERS\n"))
    vc[[i]] = vegclust(x,mobileCenters=c,nstart=nstart,...)
    if(!is.null(min.size)) {
      nsmall = sum(vc[[i]]$size< min.size)
      if(nsmall>0) {
        cat("At least one cluster was too small. Stopping.")
        if(i==1) vc<-vector("list",length=0)
        else vc<-vc[1:(i-1)]
        break
      }
    }
    i = i+1
  }
  class(vc)<-"mvegclust"
  return(vc)
}


#' @rdname hier.vegclust
#' @export
hier.vegclustdist<-function(x, hclust, cmin=2,cmax=20, min.size = NULL,  verbose=TRUE,...) {
  vc = vector("list",length=cmax-cmin +1)
  #   xs = x[row.names(x)%in% hclust$labels,] #Select those objects used for hierarchical clustering
  i=1
  for(c in cmin:cmax) {
    if(verbose) cat(paste("PROCESSING",c,"MOBILE CLUSTERS\n"))
    if(c==1) {
      memb=as.memb(rep(1, ncol(as.matrix(x))))
    }
    else memb = as.memb(cutree(hclust,k=c))
    vc[[i]] = vegclustdist(x,mobileMemb=memb,...)
    if(!is.null(min.size)) {
      nsmall = sum(vc[[i]]$size< min.size)
      if(nsmall>0) {
        cat("At least one cluster was too small. Stopping.")
        if(i==1) vc<-vector("list",length=0)
        else vc<-vc[1:(i-1)]
        break
      }
    }
    i = i+1
  }
  class(vc)<-"mvegclust"
  return(vc)
}


#' @rdname hier.vegclust
#' @export
random.vegclustdist<-function(x, cmin=2, cmax=20, nstart=10,  min.size = NULL, verbose=TRUE,...) {
  nvc = cmax-cmin +1
  vc = vector("list",nvc)
  i=1
  for(c in cmin:cmax) {
    if(verbose) cat(paste("PROCESSING",c,"MOBILE CLUSTERS\n"))
    vc[[i]] = vegclustdist(x,mobileMemb=c,nstart=nstart,...)
    if(!is.null(min.size)) {
      nsmall = sum(vc[[i]]$size< min.size)
      if(nsmall>0) {
        cat("At least one cluster was too small. Stopping.")
        if(i==1) vc<-vector("list",length=0)
        else vc<-vc[1:(i-1)]
        break
      }
    }
    i = i+1
  }
  class(vc)<-"mvegclust"
  return(vc)
}