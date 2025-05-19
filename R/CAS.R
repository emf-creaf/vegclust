#' Cumulative abundance surface (CAS)
#' 
#' Functions to calculate cumulative abundance surfaces (CASs), to build matrices from them, and to summarize several CASs.
#'
#' @param x An object of class 'doublestratifiedvegdata' (see function \code{\link{stratifyvegdata}}).
#' @param transform A function or the name of a function to be applied to each cumulative abundance value.
#' @param verbose A logical flag to indicate extra output.
#'
#' @details
#' Function \code{CAS} replaces the abundance value of each combination of size classes by the sum of abundances in this and larger size classes. This creates a cumulative abundance surface (similar to a bivariant cummulative distribution function). Function \code{CASmargin} takes an object of class '\code{CAS}' and returns an object of class '\code{CAP}' that corresponds marginal profile in either the primary or the secondary size classes. Function \code{CAS2matrix} takes an object of class '\code{CAS}' and returns a data matrix, where values differ depending on parameter \code{type}: (1) \code{type="cummulative"} simply reshapes the '\code{CAS}' object (a list) into a matrix with as many rows as plot records and where columns are organized in blocks (there are as many blocks as species and each block has as many columns as combinations of size classes); (2) \code{type="total"} returns a plot-by-species matrix where each value is the total abundance of the species in the plot (i.e. the CAS value at the ground level). When provided, \code{classWeights} are used to weight size classes of the cumulative abundance surfaces (for (1) only). Function \code{CAScenters} calculates the average abundance surface for a set of plot records. If \code{y} is a factor, it is used to speficy groups of samples for which average profiles are to be calculated. If \code{y} is an object of class '\code{\link{vegclust}}' then the function returns the CAS centroids or medoids corresponding to the clustering result. Function \code{CASquantile} calculates a quantile surface for a set of CASs. The usage of \code{y} is the same as for \code{CAScenters}.
#' 
#' @returns
#' Function \code{CAS} returns an object of class '\code{CAS}', similar to objects of class '\code{doublestratifiedvegdata}' but where abundance values of upper size classes have beed added to those of lower size classes. Function \code{CAS2matrix} returns a matrix with species as rows (columns depend on the value of \code{type}). Functions \code{CAScenters} and \code{CASquantile} return an object of class '\code{CAS}'.
#' 
#' @references
#' De \enc{Cáceres}{Caceres}, M., Legendre, P. & He, F. (2013) Dissimilarity measurements and the size structure of ecological communities. Methods in Ecology and Evolution 4: 1167-1177.
#' 
#' De \enc{Cáceres}{Caceres}, M., Coll, L., \enc{Martín-Alcón}{Martin-Alcon}, S., \enc{González-Olabarria}{Gonzalez-Olabarria}, J.R. (submitted) A general method for the classification of forest stands using structure and composition.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF.
#' 
#' @seealso \code{\link{stratifyvegdata}}, \code{\link{plot.CAS}}, \code{\link{vegdiststruct}}
#' 
#' @export
#'
#' @name CAS
#' 
#' @examples
#' ## Load tree data
#' data(treedata)
#' 
#' ## Define stratum thresholds (4 strata)
#' heights <- seq(0,4, by=0.5)
#' diameters <- seq(0,2, by=0.5)
#' 
#' ## Stratify tree data using heights and diameters as structural variables
#' X <- stratifyvegdata(treedata, sizes1=heights, sizes2=diameters, plotColumn="plotID",
#'                     speciesColumn="species", size1Column="height", size2Column="diam",
#'                     counts=TRUE)
#' X[[2]]
#' 
#' ## Build cummulative abundance surface
#' Y <- CAS(X)
#' 
#' Y[[2]]
#' 
#' ##  Extracts the first and second marginal (i.e. CAP on heights or diameters respectively)
#' Y.M1 <- CASmargin(Y, margin = 1)
#' Y.M1[[2]]
#' 
#' Y.M2 <- CASmargin(Y, margin = 2)
#' Y.M2[[2]]
#' 
#' ##  For comparison we calculate the same profiles using the stratifyvegdata and CAP functions
#' Y1 <- CAP(stratifyvegdata(treedata, sizes1=heights, plotColumn="plotID",
#'                          speciesColumn="species", size1Column="height",
#'                          counts=TRUE))                   
#' Y1[[2]]
#' Y2 <- CAP(stratifyvegdata(treedata, sizes1=diameters, plotColumn="plotID",
#'                          speciesColumn="species", size1Column="diam",
#'                          counts=TRUE))                   
#' Y2[[2]]
#' 
#' ##  Compare Y.M1[[2]] with Y1[[2]] and Y.M2[[2]] with Y2[[2]]
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

#' @rdname CAS
#' @param CAS An object of class '\code{CAS}'.
#' @param type The type of information that the resulting matrix should contain (either \code{"cummulative"} or \code{"total"}).
#' @param classWeights A numerical matrix containing the weight for each combination of size classes. If \code{NULL}, then all classes are assumed to have the same weight.
#' @export
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

#' @rdname CAS
#' @param y A vector used as a factor to calculate average or quantile surfaces per each level. Alternatively, an object of class \code{\link{vegclust}} for which CAS centroids or medoids are desired.
#' @export
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

#' @rdname CAS
#' @export
#' @param margin Indicates whether marginalization should be done in primary (\code{margin = 1}) or secondary (\code{margin = 2}) size classes.
CASmargin<-function(CAS, margin=1, verbose=FALSE) {
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

#' @rdname CAS
#' @param q Probability value for which the quantile is desired. By default the median is given.
#' @export
CASquantile<-function(CAS, q = 0.5, y = NULL) {
  quantileCAS<-function(x, q = 0.5) {
    res = x[[1]]
    v = numeric(length(x))
    d = dim(res)
    for(r in 1:d[1]) {
      for(c in 1:d[2]) {
        for(s in 1:d[3]) {
          for(i in 1:length(x)) v[i] = x[[i]][r,c,s]
          res[r,c,s] = quantile(v, probs = q)
        }
      }
    }
    return(res)
  }
  if(!is.null(y)) {
    if(is.vector(y) || is.factor(y)) {
      CC = lapply(split(CAS, as.factor(y)), FUN = quantileCAS, q)
    }    
  } else {
    CC = list(qCAS=quantileCAS(CAS, q))
  }
  if(inherits(CAS, "CAS")) class(CC)<-c("CAS","list")
  return(CC)  
}