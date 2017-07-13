CAS<-function(x, transform=NULL, verbose=FALSE) {
  cas<-function(x, verbose=FALSE) {
    if(verbose) cat(".")
    a = as.array(x)
    nsp = dim(a)[1]
    n2 = dim(a)[2]
    n3 = dim(a)[3]
    for(sp in 1:nsp) {
      if(n2>1 || n3>1) {
        for(i in n2:1) {
          for(j in n3:1) {
            a[sp,i,j] = sum(x[sp,i:n2, j:n3])
          }
        }
      }      
    }
    return(a)
  }
  if(!inherits(x,"doublestratifiedvegdata")) stop("Input should be of class 'doublestratifiedvegdata'")
  Y = lapply(x, FUN=cas, verbose=verbose)
  if(!is.null(transform)) Y = lapply(Y, FUN=transform)
  names(Y)<-names(x)  
  class(Y)<-c("CAS","list")
  return(Y)
}