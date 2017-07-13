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