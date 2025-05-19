#' Cumulative abundance profile (CAP)
#'
#' Functions to calculate cumulative abundance profiles (CAPs), to build matrices from them, and to summarize several profiles.
#' 
#' @param x A stratified vegetation data set (see function \code{\link{stratifyvegdata}}).
#' @param transform A function or the name of a function to be applied to each cumulative abundance value.
#' @param verbose A logical flag to indicate extra output.
#'
#' @returns
#' Function \code{CAP} returns an object of class '\code{CAP}', similar to objects of class '\code{stratifiedvegdata}' but where abundance values of upper size classes have beed added to those of lower size classes. Function \code{CAP2matrix} returns a matrix with species as rows (columns depend on the value of \code{type}). Functions \code{CAPcenters} and \code{CAPquantile} return an object of class '\code{CAP}'.
#' 
#' @details
#' Function \code{CAP} replaces the abundance value of a size class by the sum of abundances in this and larger size classes (strata). Thus, upper classes contain smaller abundance values than lower classes, creating a cumulative abundance profile. Function \code{CAP2matrix} takes an object of class '\code{CAP}' and returns a data matrix, where values differ depending on parameter \code{type}: (1) \code{type="cumulative"} simply reshapes the '\code{CAP}' object (a list) into a matrix with as many rows as plot records and where columns are organized in blocks (there are as many blocks as species and each block has as many columns as size classes); (2) \code{type="total"} returns a plot-by-species matrix where each value is the total abundance of the species in the plot (i.e. the CAP value at the ground level); (3) \code{type="volume"} returns a plot-by-species matrix where each value is the sum of CAP values across size classes (a measure of the "volume" occupied by the species in the plot). When provided, \code{classWeights} are used to weight size classes of the cumulative abundance profiles (for (1) and (3) only). Function \code{CAPcenters} calculates the average abundance profile for a set of plot records. If \code{y} is a factor, it is used to speficy groups of samples for which average profiles are to be calculated. If \code{y} is an object of class '\code{\link{vegclust}}' then the function returns the CAP centroids or medoids corresponding to the clustering result. Function \code{CAPquantile} calculates a quantile profile for a set of CAPs. The usage of \code{y} is the same as for \code{CAPcenters}.
#' 
#' @references
#' De \enc{Cáceres}{Caceres}, M., Legendre, P. & He, F. (2013) Dissimilarity measurements and the size structure of ecological communities. Methods in Ecology and Evolution 4: 1167-1177.
#' 
#' De \enc{Cáceres}{Caceres}, M., Coll, L., \enc{Martín-Alcón}{Martin-Alcon}, S., \enc{González-Olabarria}{Gonzalez-Olabarria}, J.R. (submitted) A general method for the classification of forest stands using structure and composition.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF.
#' 
#' @seealso \code{\link{stratifyvegdata}}, \code{\link{plot.CAP}}, \code{\link{vegdiststruct}}
#' 
#' @export
#'
#' @name CAP
#' @examples
#' ## Load stratified data
#' data(medreg)
#' 
#' ## Check that 'medreg' has correct class
#' class(medreg)
#' 
#' ## Look at the data for the third plot
#' medreg[[3]]
#' 
#' ## Create cumulative abundance profile (CAP) for each plot
#' medreg.CAP <- CAP(medreg)
#' 
#' ## Look at the profile of the third plot
#' medreg.CAP[[3]]
#' 
#' ## Create matrix with species abundances
#' medreg.X <- CAP2matrix(medreg.CAP, type="total")
#' head(medreg.X)
#' 
#' ## Generate and plot average profile
#' average.CAP <- CAPcenters(medreg.CAP)
#' plot(average.CAP)
#' 
#' ## Generate and plot median profile
#' median.CAP <- CAPquantile(medreg.CAP, q = 0.5)
#' plot(median.CAP)
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

#' @rdname CAP
#' @param CAP An object of class '\code{CAP}'.
#' @param type The type of information that the resulting matrix should contain. Either \code{"profile"}, \code{"abundance"} or \code{"volume"}.
#' @param classWeights A numerical vector containing the weight for size class. If \code{NULL}, then all classes are assumed to have the same weight.
#' @export
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
  if(is.null(spnames)) spnames = 1:nspecies
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

#' @rdname CAP
#' @export
#' @param y A vector used as a factor to calculate average or quantile profiles per each level. Alternatively, an object of class \code{\link{vegclust}} for which CAP centroids or medoids are desired.
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

#' @rdname CAP
#' @param q Probability value for which the quantile is desired. By default the median is given.
#' @export
CAPquantile<-function(CAP, q = 0.5, y = NULL) {
  quantileCAP<-function(x, q = 0.5) {
    res = x[[1]]
    v = numeric(length(x))
    for(r in 1:nrow(res)) {
      for(c in 1:ncol(res)) {
        for(i in 1:length(x)) v[i] = x[[i]][r,c]
        res[r,c] = quantile(v, probs = q)
      }
    }
    return(res)
  }
  if(!is.null(y)) {
    if(is.vector(y) || is.factor(y)) {
      CC = lapply(split(CAP, as.factor(y)), FUN = quantileCAP, q)
    }    
  } else {
    CC = list(qCAP=quantileCAP(CAP, q))
  }
  if(inherits(CAP, "CAP")) class(CC)<-c("CAP","list")
  else if(inherits(CAP, "stratifiedvegdata")) class(CC)<-c("stratifiedvegdata","list")
  return(CC)  
}