#' Structural and compositional dissimilarity
#' 
#' Function to calculate the dissimilarity between ecological communities taking into account both their composition and the size of organisms.
#'
#' @param x A stratified vegetation data set (see function \code{\link{stratifyvegdata}}), a set of cummulative abundance profiles (see function \code{\link{CAP}}) or a set of cummulative abundance surfaces (see function \code{\link{CAS}}).
#' @param y A second stratified vegetation data set (see function \code{\link{stratifyvegdata}}), a second set of cummulative abundance profiles (see function \code{\link{CAP}}) or a second set of cummulative abundance surfaces (see function \code{\link{CAS}}) against which object \code{x} should be compared.
#' @param paired Only relevant when \code{y != NULL}. If \code{paired = TRUE} pairwise comparisons are calculated between elements in \code{x} and \code{y} (and \code{x} and \code{y} need to be of the same length). If \code{paired = FALSE} then all objects in \code{x} are compared to all objects in \code{y}.
#' @param type Whether dissimilarities between pairs of sites should be calculated from differences in cummulative abundance (\code{"cumulative"}), in total abundance (\code{"total"}) or in volumes of cumulative abundance profiles (\code{"volume"}).
#' @param method The dissimilarity coefficient to calculate (see details).
#' @param transform A function or the name of a function to be applied to each cumulative abundance value.
#' @param classWeights A numerical vector or a matrix containing the weight of each size class or combination of size classes (see functions \code{\link{CAP2matrix}} and \code{\link{CAS2matrix}}). If \code{NULL}, then the function assumes classes of equal weight.
#'
#' @details
#' The six different coefficients available are described in De Caceres et al. (2013): (1) \code{method="bray"} for percentage difference (alias Bray-Curtis dissimilarity); (2) \code{method="ruzicka"} for Ruzicka index (a generalization of Jaccard); (3) \code{method="kulczynski"} for the Kulczynski dissimilarity index; (4) \code{method="ochiai"} for the complement of a quantitative generalization of Ochiai index of similarity; (5) \code{method="canberra"} for the Canberra index (Adkins form); (6) \code{method="relman"} for the relativized Manhattan coefficient (Whittaker's index of association). Currently, the function also supports (7) \code{method="manhattan"} for the city block metric.
#' 
#' @returns
#' Returns an object of class '\code{\link{dist}}'.
#' 
#' @seealso \code{\link{stratifyvegdata}}, \code{\link[vegan]{vegdist}}
#' 
#' @references
#' De \enc{Cáceres}{Caceres}, M., Legendre, P. & He, F. (2013) Dissimilarity measurements and the size structure of ecological communities. Methods in Ecology and Evolution 4: 1167-1177.
#' 
#' @export
#'
#' @examples
#' ## Load stratified data
#' data(medreg)
#' 
#' ## Check that 'medreg' has correct class
#' class(medreg)
#' 
#' ## Create cumulative abundance profile (CAP) for each plot
#' medreg.CAP <- CAP(medreg)
#' 
#' ## Create dissimilarity (percentage difference) matrix using profiles
#' medreg.D <- vegdiststruct(medreg, method="bray")
#' 
#' ## Create dissimilarity (percentage difference) matrix using abundances
#' medreg.D2 <- vegdiststruct(medreg, method="bray", type="total")
#' 
#' ## Calculate correlation
#' cor(as.vector(medreg.D), as.vector(medreg.D2))
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