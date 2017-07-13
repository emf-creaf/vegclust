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