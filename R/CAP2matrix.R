CAP2matrix<-function(CAP, type="cumulative", classWeights=NULL) {
  type= match.arg(type,c("cumulative","total","volume"))  
  n=length(CAP)
  if(!inherits(CAP,"CAP")) stop("Wrong input data. Use function CAP() first.")
  if(!is.null(classWeights)) {
    if(ncol(CAP[[1]])!=length(classWeights)) stop("Number of size classes in stratified species data does not match the number of elements in the vector of class weights")
  }
  else classWeights = rep(1,ncol(CAP[[1]])) #All classes have equal width
  nstrata = length(classWeights)
  nspecies = nrow(CAP[[1]])
  spnames = row.names(CAP[[1]])
  if(type=="cumulative") {
    m = matrix(0, nrow=n,ncol=nstrata*nspecies)
    for(i in 1:n) {
      m[i,] = as.vector(t(CAP[[i]])*classWeights)
    }    
    colnames(m) = paste(as.character(gl(nspecies,nstrata, labels=spnames)),1:nstrata, sep="_")
  } else if(type=="volume") {
    m = matrix(0, nrow=n,ncol=nspecies)
    for(i in 1:n) {
      wm = t(CAP[[i]])*classWeights
      m[i,] = colSums(wm)
    }
    colnames(m) = spnames
  } else if(type=="total") {
    m = matrix(0, nrow=n,ncol=nspecies)
    for(i in 1:n) {
      for(j in 1:nspecies){
        m[i,j] = CAP[[i]][j,1]
      }
    }    
    colnames(m) = spnames
  }
  rownames(m) = names(CAP)
  return(m)
}