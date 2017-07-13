CAScenters<-function(CAS, y=NULL) {
  averageCAS<-function(x) {
    avccf = x[[1]]
    if(length(x)>1) {
      for(i in 2:length(x)) {
        avccf = avccf + x[[i]]
      }
      avccf = avccf/length(x)
    }
    return(avccf)
  }
  sumCAS<-function(x) {
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
      CC = lapply(split(CAS, as.factor(y)), FUN = averageCAS)
    } else if(inherits(y,"vegclust")) {
      if(y$method %in% c("KMdd","FCMdd","HNCdd","NCdd","PCMdd")) {
        mi = c(y$mobileCenters, y$fixedCenters)
        CC = CAS[mi]
      } else {
        memb = y$memb
        CC = vector("list", ncol(memb))
        if(nrow(memb)!=length(CAS)) stop("The number of plots in CAS has to be equal to the classified elements in y")
        for(cl in 1:ncol(memb)) {
          CASmemb = CAS
          s = sum(memb[,cl])
          for(i in 1:nrow(memb)) CASmemb[[i]] = CASmemb[[i]]*memb[i,cl]
          CC[[cl]]<-sumCAS(CASmemb)/s
          names(CC)[cl]<-colnames(memb)[cl]
        }      
      }
    }
  } else {
    CC = list(meanCAS=averageCAS(CAS))
  }
  class(CC)<-c("CAS","list")
  return(CC)  
}