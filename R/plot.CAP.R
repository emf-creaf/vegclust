#' Draws cummulative abundance profiles
#'
#' Create plots used to inspect one or more cumulative abundance profiles.
#' 
#' @param x An object returned from function \code{\link{CAP}} or an object of class \code{stratifiedvegdata} (see documentation for function \code{\link{stratifyvegdata}}).
#' @param sizes A vector containing the size values associated to each size class. If \code{NULL} the y-axis will be defined using the size class order in \code{x}.
#' @param species A vector of strings indicating the species whose profile is to be drawn. If \code{NULL} all species are plotted.
#' @param plots A vector indicating the plot records whose profile is to be drawn. Can be a \code{character} vector (for plot names), a \code{numeric} vector (for plot indices) or a \code{logical} vector (for TRUE/FALSE selection). If \code{NULL} all plot records are plotted.
#' @param switchAxes A flag indicating whether ordinate and abscissa axes should be interchanged.
#' @param add A flag indicating whether profiles should be drawn on top of current drawing area. If \code{add=FALSE} a new plot is created.
#' @param drawAxes A flag indicating whether axes should be drawn.
#' @param xlab String label for the x axis.
#' @param ylab String label for the y axis.
#' @param type Type of plot to be drawn ("p" for points, "l" for lines, "s" for steps, ...).
#' @param ... Additional plotting parameters.
#'
#' @references
#' De \enc{Cáceres}{Caceres}, M., Legendre, P. & He, F. (2013) Dissimilarity measurements and the size structure of ecological communities. Methods in Ecology and Evolution 4: 1167-1177.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' 
#' @seealso \code{\link{CAP}}
#' 
#' @name plot.CAP
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
#' ## Draw the stratified data and profile corresponding to the third plot
#' plot(medreg, plots="3")
#' plot(medreg.CAP, plots="3")
#' 
#' ## Look at the plot and CAP of the same plot
#' medreg[["3"]]
#' medreg.CAP[["3"]]
plot.CAP<-function(x, sizes=NULL, species=NULL, plots=NULL, switchAxes=FALSE, add=FALSE, drawAxes=TRUE, xlab="", ylab="", type="s",...) {
  if(!is.null(plots)) x = x[plots]
  if(length(x)==0) stop("No site samples to draw. Check your selection.")
  spnames = row.names(as.data.frame(x[[1]]))
  if(length(spnames)==0) stop("No species to draw. Check your selection.")
  if(!is.null(species)) spnames = spnames[spnames %in% species]
  if(length(spnames)==0) stop("No species to draw. Check your selection.")
  if(is.null(sizes)) sizes = 1:ncol(x[[1]])
  maxval = 0
  for(i in 1:length(x)) {
    maxval = max(maxval,max(as.matrix(x[[i]])))
  }
  if(!add) {
    if(switchAxes) {
      plot(x=as.numeric(c(0,rep(maxval, length(sizes)))), y=c(0,sizes), 
                      type="n", axes=FALSE,xlab = xlab, ylab=ylab,...)
    }
    else {
      plot(x=c(0,sizes), y=as.numeric(c(0,rep(maxval, length(sizes)))),
              type="n", axes=FALSE,xlab = xlab, ylab=ylab,...)
    }
    if(drawAxes) {
      axis(1,pos=0)
      axis(2,pos=0)
    }
  }
  for(i in 1:length(x)) {
    cap = x[[i]]
    cap  = cap[rownames(cap) %in% spnames,, drop=FALSE]
    cap = cbind(cap,0)
    if(switchAxes) matlines(x=t(cap), y=c(0,sizes), type=type,...)
    else matlines(x=c(0,sizes), y=t(cap), type=type,...)
  }
}

#' @rdname plot.CAP
#' @export
plot.stratifiedvegdata<-function(x, sizes=NULL, species=NULL, plots=NULL, switchAxes=FALSE, add=FALSE, drawAxes=TRUE, xlab="", ylab="", type="s",...) {
  if(!is.null(plots)) x = x[plots]
  if(length(x)==0) stop("No site samples to draw. Check your selection.")
  spnames = row.names(as.data.frame(x[[1]]))
  if(length(spnames)==0) stop("No species to draw. Check your selection.")
  if(!is.null(species)) spnames = spnames[spnames %in% species]
  if(length(spnames)==0) stop("No species to draw. Check your selection.")
  if(is.null(sizes)) sizes = 1:ncol(x[[1]])
  maxval = 0
  for(i in 1:length(x)) {
    maxval = max(maxval,max(as.matrix(x[[i]])))
  }
  if(!add) {
    if(switchAxes) {
      plot(x=as.numeric(c(0,rep(maxval, length(sizes)))), y=c(0,sizes), 
           type="n", axes=FALSE,xlab = xlab, ylab=ylab,...)
    }
    else {
      plot(x=c(0,sizes), y=as.numeric(c(0,rep(maxval, length(sizes)))),
           type="n", axes=FALSE,xlab = xlab, ylab=ylab,...)
    }
    if(drawAxes) {
      axis(1,pos=0)
      axis(2,pos=0)
    }
  }
  for(i in 1:length(x)) {
    cap = x[[i]]
    cap  = cap[rownames(cap) %in% spnames,, drop=FALSE]
    cap = cbind(cap,0)
    if(switchAxes) matlines(x=t(cap), y=c(0,sizes), type=type,...)
    else matlines(x=c(0,sizes), y=t(cap), type=type,...)
  }
}