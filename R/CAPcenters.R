CAPcenters<-function(CAP, y=NULL) {
  averageCAP<-function(x) {
    avccf = x[[1]]
    if(length(x)>1) {
      for(i in 2:length(x)) {
        avccf = avccf + x[[i]]
      }
      avccf = avccf/length(x)
    }
    return(avccf)
  }
  sumCAP<-function(x) {
    sccf = x[[1]]
    if(length(x)>1) {
      for(i in 2:length(x)) {
        sccf = sccf + x[[i]]
      }
    }
    return(sccf)
  }  
  if(!is.null(y)) {
    if(is.vector(y) || is.factor(y)) {
      CC = lapply(split(CAP, as.factor(y)), FUN = averageCAP)
    } else if(inherits(y,"vegclust")) {
      if(y$method %in% c("KMdd","FCMdd","HNCdd","NCdd","PCMdd")) {
        mi = c(y$mobileCenters, y$fixedCenters)
        CC = CAP[mi]
      } else {
        memb = y$memb
        CC = vector("list", ncol(memb))
        if(nrow(memb)!=length(CAP)) stop("The number of plots in CAP has to be equal to the classified elements in y")
        for(cl in 1:ncol(memb)) {
          CAPmemb = CAP
          s = sum(memb[,cl])
          for(i in 1:nrow(memb)) CAPmemb[[i]] = CAPmemb[[i]]*memb[i,cl]
          CC[[cl]]<-sumCAP(CAPmemb)/s
        }
        names(CC)<-colnames(y$dist2clusters)
      }
    }    
  } else {
    CC = list(meanCAP=averageCAP(CAP))
  }
  if(inherits(CAP, "CAP")) class(CC)<-c("CAP","list")
  else if(inherits(CAP, "stratifiedvegdata")) class(CC)<-c("stratifiedvegdata","list")
  return(CC)  
}