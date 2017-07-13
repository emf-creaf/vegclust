vegdiststruct<-function(x, y=NULL, paired=FALSE, type="cumulative", method="bray", transform=NULL, classWeights=NULL) {
  method= match.arg(method,c("bray","ruzicka","kulczynski","ochiai", "canberra","relman","manhattan"))
  type= match.arg(type,c("cumulative","total","volume"))
  if(inherits(x,"stratifiedvegdata")) {
    x = CAP(x, transform=transform)
  } else if(inherits(x,"doublestratifiedvegdata")) {
    x = CAS(x, transform=transform)
  } else if(inherits(x,"CAP")) {
    if(!is.null(transform)) x = lapply(x, FUN=transform)
    class(x)<-c("list","CAP")
  } else if(inherits(x,"CAS")) {
    if(!is.null(transform)) x = lapply(x, FUN=transform)
    class(x)<-c("list","CAS")
  } else{
    stop("Wrong data type for 'x'. Please use a 'stratifiedvegdata' object, a 'doublestratifiedvegdata' object, a CAP' object or a 'CAS' object.")
  }
  if(!is.null(y)) {
    if(inherits(y,"stratifiedvegdata")) {
      y = CAP(y, transform=transform)
    } else if(inherits(y,"doublestratifiedvegdata")) {
      y = CAS(y, transform=transform)
    } else if(inherits(y,"CAP")) {
      if(!is.null(transform)) y = lapply(y, FUN=transform)
      class(y)<-c("list","CAP")
    } else if(inherits(y,"CAS")) {
      if(!is.null(transform)) y = lapply(y, FUN=transform)
      class(y)<-c("list","CAS")
    } else{
      stop("Wrong data type for 'y'. Please use a 'stratifiedvegdata' object, a 'doublestratifiedvegdata' object, a CAP' object or a 'CAS' object.")
    }    
    if(length(y)!=length(x) && paired) stop("The number of sites in 'x' and 'y' must be the same.")
  } 
  # Compare all plots in 'x' between them (squared symmetric dissimilarity matrix)
  if(is.null(y)) {
    if(inherits(x, "CAP")) Y = CAP2matrix(x,type=type, classWeights=classWeights)
    else if(inherits(x, "CAS")) Y = CAS2matrix(x,type=type, classWeights=classWeights)
    n = nrow(Y)
    res = matrix(0,n,n)
    Ysums = rowSums(Y)
    if(method=="bray") {
      for(i in 2:n) {
        for(j in 1:(n-1)) {
          A = sum(pmin(Y[i,],Y[j,]))
          res[i,j] = 1-(2*A/(Ysums[i]+Ysums[j]))
        }
      }
    } else if(method=="relman") {
      return(dist(decostand(Y,method="total"), method="manhattan"))
    } else if(method=="manhattan") {
      return(dist(Y, method="manhattan"))
    } else if(method=="ruzicka") {
      for(i in 2:n) {
        for(j in 1:(n-1)) {
          A = sum(pmin(Y[i,],Y[j,]))
          res[i,j] = 1-(A/(Ysums[i]+Ysums[j]-A))
        }
      }
    } else if(method=="kulczynski") {
      for(i in 2:n) {
        for(j in 1:(n-1)) {
          A = sum(pmin(Y[i,],Y[j,]))
          res[i,j] = 1-0.5*((A/Ysums[i])+(A/Ysums[j]))
        }
      }
    } else if(method=="ochiai") {
      for(i in 2:n) {
        for(j in 1:(n-1)) {
          A = sum(pmin(Y[i,],Y[j,]))
          res[i,j] = 1-(A/sqrt(Ysums[i]*Ysums[j]))
        }
      }
    } else if(method=="canberra") {
      p = nrow(x[[1]]) #Number of species
      s = ncol(x[[1]]) #Number of strata
      for(i in 2:n) {
        for(j in 1:(n-1)) {
          num = colSums(matrix((abs(Y[i,]-Y[j,])), nrow=s, ncol=p))
          den = colSums(matrix((Y[i,]+Y[j,]), nrow=s, ncol=p))
          sel = (den>0)
          pp = sum(sel)
          num = num[sel]
          den = den[sel]
          res[i,j] = sum(num/den)/pp
        }
      }
    }
    rownames(res) = colnames(res) = names(x)
    return(as.dist(res))    
  } else if(!is.null(y) && paired==TRUE) {
    if(inherits(x,"CAP")) X = CAP2matrix(x,type=type, classWeights=classWeights)
    else if(inherits(x, "CAS")) X = CAS2matrix(x,type=type, classWeights=classWeights)
    n = nrow(X)
    if(inherits(y,"CAP")) Y = CAP2matrix(y,type=type, classWeights=classWeights)
    else if(inherits(y,"CAS")) Y = CAS2matrix(y,type=type, classWeights=classWeights)
    res = rep(0,n)
    Xsums = rowSums(X)
    Ysums = rowSums(Y)
    if(method=="bray") {
      for(i in 1:n) {
        res[i] = 1-(2*sum(pmin(X[i,],Y[i,]))/(Xsums[i]+Ysums[i]))
      }
    } else if(method=="relman") {
      for(i in n) {
        res[i] = sum(abs((X[i,]/Xsums[i])-(Y[i,]/Ysums[i])))
      }
    } else if(method=="manhattan") {
      for(i in n) {
        res[i] = sum(abs(X[i,]-Y[i,]))
      }
    } else if(method=="ruzicka") {
      for(i in 1:n) {
        A = sum(pmin(X[i,],Y[i,]))
        res[i] = 1-(A/(Xsums[i]+Ysums[i]-A))
      }
    } else if(method=="kulczynski") {
      for(i in 1:n) {
        A = sum(pmin(X[i,],Y[i,]))
        res[i] = 1-0.5*((A/Ysums[i])+(A/Ysums[i]))
      }
    } else if(method=="ochiai") {
      for(i in 1:n) {
        res[i] = 1-(sum(pmin(X[i,],Y[i,]))/sqrt(Ysums[i]*Ysums[i]))
      }
    } else if(method=="canberra") {
      p = nrow(x[[1]]) #Number of species
      s = ncol(x[[1]]) #Number of strata
      for(i in 1:n) {
        num = colSums(matrix((abs(X[i,]-Y[i,])), nrow=s, ncol=p))
        den = colSums(matrix((X[i,]+Y[i,]), nrow=s, ncol=p))
        sel = (den>0)
        pp = sum(sel)
        num = num[sel]
        den = den[sel]
        res[i] = sum(num/den)/pp
      }
    }
    return(res)
    
  # Rectangular matrix  
  } else if(!is.null(y) && paired==FALSE) {
    if(inherits(x,"CAP")) X = CAP2matrix(x,type=type, classWeights=classWeights)
    else if(inherits(x, "CAS")) X = CAS2matrix(x,type=type, classWeights=classWeights)
    nX = nrow(X)
    if(inherits(y,"CAP")) Y = CAP2matrix(y,type=type, classWeights=classWeights)
    else if(inherits(y,"CAS")) Y = CAS2matrix(y,type=type, classWeights=classWeights)
    nY = nrow(Y)
    res = matrix(0,nX,nY)
    rownames(res)<-rownames(X)
    colnames(res)<-rownames(Y)
    Xsums = rowSums(X)
    Ysums = rowSums(Y)
    if(method=="bray") {
      for(i in 1:nX) {
        for(j in 1:nY) {
          res[i,j] = 1-(2*sum(pmin(X[i,],Y[j,]))/(Xsums[i]+Ysums[j]))
        }
      }
    } else if(method=="relman") {
      for(i in 1:nX) {
        for(j in 1:nY) {
          res[i,j] = sum(abs((X[i,]/Xsums[i])-(Y[j,]/Ysums[j])))
        }
      }
    } else if(method=="manhattan") {
      for(i in 1:nX) {
        for(j in 1:nY) {
          res[i,j] = sum(abs(X[i,]-Y[j,]))
        }
      }
    } else if(method=="ruzicka") {
      for(i in 1:nX) {
        for(j in 1:nY) {
          A = sum(pmin(X[i,],Y[j,]))
          res[i,j] = 1-(A/(Xsums[i]+Ysums[j]-A))
        }
      }
    } else if(method=="kulczynski") {
      for(i in 1:nX) {
        for(j in 1:nY) {
          A = sum(pmin(X[i,],Y[j,]))
          res[i,j] = 1-0.5*((A/Ysums[i])+(A/Ysums[j]))
        }
      }
    } else if(method=="ochiai") {
      for(i in 1:nX) {
        for(j in 1:nY) {
          res[i,j] = 1-(sum(pmin(X[i,],Y[j,]))/sqrt(Ysums[i]*Ysums[j]))
        }
      }
    } else if(method=="canberra") {
      p = nrow(x[[1]]) #Number of species
      s = ncol(x[[1]]) #Number of strata
      for(i in 1:nX) {
        for(j in 1:nY) {
          num = colSums(matrix((abs(X[i,]-Y[j,])), nrow=s, ncol=p))
          den = colSums(matrix((X[i,]+Y[j,]), nrow=s, ncol=p))
          sel = (den>0)
          pp = sum(sel)
          num = num[sel]
          den = den[sel]
          res[i,j] = sum(num/den)/pp
        }
      }
    }
    return(res)
  } 
}