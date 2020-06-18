CAP<-function(x, transform=NULL, verbose=FALSE) {
  cap<-function(x, verbose=FALSE) {
    if(verbose) cat(".")
    y = as.matrix(x)
    if(ncol(y)>1) {
      for(i in (ncol(y)-1):1) {
        y[,i] = y[,i]+y[,i+1]
      }
    }
    dimnames(y) <- dimnames(x)
    return(y)
  }
  if(!inherits(x,"stratifiedvegdata")) stop("Input should be of class 'stratifiedvegdata'")
  Y = lapply(x, FUN=cap, verbose=verbose)
  if(!is.null(transform)) Y = lapply(Y, FUN=transform)
  names(Y)<-names(x)
  class(Y)<-c("CAP","list")
  return(Y)
}