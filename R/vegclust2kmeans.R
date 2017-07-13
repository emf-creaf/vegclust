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