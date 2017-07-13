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