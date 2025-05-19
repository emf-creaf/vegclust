#' Concordance between two classifications
#' 
#' Computes an index to compare two classifications.
#'
#' @param x,y Classification vector or membership matrix. Alternatively, objects of type \code{\link{vegclust}} or \code{\link{vegclass}}. 
#' @param method A string vector to indicate the desired indices (see details). 
#' @param ... Additional parameters for function \code{\link{defuzzify}}, which will be called if \code{x} or \code{y} are of type \code{matrix}, \code{\link{vegclust}} or \code{\link{vegclass}}.
#'
#' @details
#' Several indices for comparison of partitions are available:
#' \itemize{
#'   \item{\code{method="Rand"}: Rand (1971) index.}
#'   \item{\code{method="adjustedRand"}: Rand index adjusted for random effects (Hubert & Arabie 1985).}
#'   \item{\code{method="Wallace"}: Wallace (1983) index (for asymmetrical comparisons). This index (and its adjusted version) is useful to quantify how much \code{x} is nested into \code{y}.}
#'   \item{\code{method="adjustedWallace"}: Wallace index adjusted for random effects (Pinto et al. 2008).}
#' }
#' 
#' @references
#' Hubert, L. & Arabie, P. (1985). Comparing partitions. Journal of Classification, 2, 193–218.
#' 
#' Pinto, F.R., Melo-Cristino, J. & Ramirez, M. (2008). A confidence interval for the wallace coefficient of concordance and its application to microbial typing methods. PLoS ONE, 3.
#' 
#' Rand, W.M. (1971). Objective Criteria for the Evaluation of Clustering Methods. Journal of the American Statistical Association, 66, 846–850.
#' 
#' Wallace, D.L. (1983). A method for comparing two hierarchical clusterings: Comment. Journal of the American Statistical Association, 78, 569–576.
#' 
#' @returns 
#' A numeric vector with the desired index values.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' 
#' @seealso \code{\link{vegclust}}, \code{\link{vegclass}}, \code{\link{defuzzify}}
#' @export
#' 
concordance<-function(x,y, method="adjustedRand",...){
  method = match.arg(method,choices=c("Rand","adjustedRand","Wallace","adjustedWallace"), 
                      several.ok=TRUE)
  if(inherits(x, "vegclust") || inherits(x, "vegclass") || inherits(x, "matrix")) x<-defuzzify(x,...)$cluster
  if(inherits(y, "vegclust") || inherits(y, "vegclass") || inherits(y, "matrix")) y<-defuzzify(y,...)$cluster
  if (length(x) != length(y)) 
    stop("arguments must be classifications of the same number of objects")
  tab <- table(x, y)
  N<-length(x)
  if (all(dim(tab) == c(1, 1))) return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  I<-numeric(length(method))
  for(i in 1:length(method)) {
    if(method[i]=="Rand") {
      I[i]<-(a+d)/(a+b+c+d)
    } else if(method[i]=="adjustedRand"){
      I[i] <- (a - (a + b) * (a + c)/(a + b + c + d))/
        ((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
    } else if(method[i]=="Wallace") {
      I[i] <- a/(a+b)
    } else if(method[i]=="adjustedWallace"){
      W <-a/(a+b)
      nj = colSums(tab)
      SID = (sum(nj*(nj-1))/(N*(N-1)))
      I[i] = (W-SID)/(1-SID)
    }
  }
  names(I)<-method
  return(I)
}