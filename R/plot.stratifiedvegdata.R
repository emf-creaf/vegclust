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