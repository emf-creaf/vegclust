CAS2matrix<-function(CAS, type="cumulative", classWeights=NULL) {
  if(!inherits(CAS,"CAS")) stop("Wrong input data. Use function CAS() first.")
  type= match.arg(type,c("cumulative","total"))  
  n=length(CAS)
  if(n<2) stop("Wrong number of plot records (should be larger than one)")
  nspecies = dim(CAS[[1]])[1]
  spnames = dimnames(CAS[[1]])[[1]]
  nsizes1 = dim(CAS[[1]])[2]
  nsizes2 = dim(CAS[[1]])[3]  
  if(!is.null(classWeights)) {
    if(nsizes1!=nrow(classWeights)) stop("Number of size1 classes in 'CAS' does not match the number of rows in matrix 'classWeights'")
    if(nsizes2!=ncol(classWeights)) stop("Number of size2 classes in 'CAS' does not match the number of columns in matrix 'classWeights'")
  }
  else {
    classWeights = matrix(1, nrow=nsizes1, ncol=nsizes2)
  } 
  if(type=="cumulative") {
    m = matrix(0, nrow=n,ncol=nsizes1*nsizes2*nspecies)
    for(i in 1:n) {
      m[i,] = as.vector(aperm(CAS[[i]],c(3,2,1)))*as.vector(classWeights)
    }    
#     colnames(m) = paste(as.character(gl(nspecies,nsizes1*nsizes2, labels=spnames)),1:nstrata, sep="_")
  } else if(type=="total") {
    m = matrix(0, nrow=n,ncol=nspecies)
    for(i in 1:n) {
      for(j in 1:nspecies){
        m[i,j] = CAS[[i]][j,1,1]
      }
    }    
    colnames(m) = spnames
  }
  rownames(m) = names(CAS)
  return(m)
}