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