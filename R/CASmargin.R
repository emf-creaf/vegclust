CASmargin<-function(CAS,margin=1, verbose=FALSE) {
  if(!inherits(CAS,"CAS")) stop("Input should be of class 'CAS'")
  margin= match.arg(as.character(margin),choices=c("1","2"))
  cas2cap<-function(x, margin="1", verbose=FALSE) {
    if(verbose) cat(".")
    if(margin=="1") {
      y = matrix(x[,,1],nrow=dim(x)[1], ncol=dim(x)[2])
      dn <-dimnames(x)
      dimnames(y)<-list(dn[[1]],dn[[2]]) 
    } else if(margin=="2"){
      y = matrix(x[,1,],nrow=dim(x)[1], ncol=dim(x)[3])
      dn <-dimnames(x)      
      dimnames(y)<-list(dn[[1]],dn[[3]]) 
    }
    return(y)
  }
  Y = lapply(CAS, FUN=cas2cap, margin=margin, verbose=verbose)
  names(Y)<-names(CAS)  
  class(Y)<-c("CAP","list")
  return(Y)
}