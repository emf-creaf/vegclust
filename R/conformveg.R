#' Conform two community data tables
#' 
#' Conforms two community data tables to have the same set of columns (species)
#'
#' @param x Community data, a site-by-species matrix.
#' @param y Community data, a site-by-species matrix.
#' @param fillvalue The value to be used to fill new entries in inflated matrices.
#' @param verbose Displays information about the number of species shared between \code{x} and \code{y}, as well as the number of species that are in one of the data tables but not in the other.
#'
#' @details
#' This function adds to \code{x} as many new columns as columns of \code{y} that are not in \code{x}. The same is done for \code{y}, so the two tables have the same set of columns when they are returned. 
#' 
#' @returns
#' A list with the two inflated matrices \code{x} and \code{y}.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF.
#' 
#' @seealso \code{\link{vegclust}}, \code{\link{vegclass}}
#' 
#' @export
#'
#' @examples
#' ## Loads data (38 columns and 33 species)
#' data(wetland)
#' dim(wetland)
#' 
#' ## Splits wetland data into two matrices of 30x27 and 11x22
#' wetland.30 <- wetland[1:30,]
#' wetland.30 <- wetland.30[,colSums(wetland.30)>0]
#' dim(wetland.30)
#' wetland.11 <- wetland[31:41,]
#' wetland.11 <- wetland.11[,colSums(wetland.11)>0] 
#' dim(wetland.11)
#' 
#' ## Conforms the two matrices so they can eventually be merged
#' wetland.cf <- conformveg(wetland.30, wetland.11)
#' dim(wetland.cf$x)
#' dim(wetland.cf$y)
#' names(wetland.cf$x)==names(wetland.cf$y)
conformveg<-function(x,y, fillvalue=0, verbose=FALSE) {
	nx = names(x)
	ny = names(y)
	if(is.null(nx) || is.null(ny)) {
		stop("x and y must have both column names")
	}
	nn = ny[!(ny %in% nx)]
	nm = c(nx,nn)
	if(verbose) {
		cat(paste(sum(ny %in% nx)," names in Y that are in X."))
		cat(paste(sum(!(ny %in% nx))," names in Y that are NOT in X: "))
		nynx = ny[!(ny %in% nx)]
		for(i in 1:length(nynx)) cat(paste(nynx[i],"-"))
		cat("\n")
		cat(paste(sum(nx %in% ny)," names in X that are in Y"))
		cat(paste(sum(!(nx %in% ny))," names in X that are NOT in Y: "))
		nxny = nx[!(nx %in% ny)]
		for(i in 1:length(nxny)) cat(paste(nxny[i],"-"))
		cat("\n")
	}
	
	xinf = data.frame(matrix(fillvalue,nrow=nrow(x),ncol=length(nm)))
	xinf[,1:length(nx)] = x	
	row.names(xinf) = row.names(x)
	names(xinf) = nm
	
	yinf = data.frame(matrix(fillvalue,nrow=nrow(y),ncol=length(nm)))
	#Shared columns
	sel = which(ny %in% nx)
	for(i in 1:length(sel)) yinf[,which(nx %in% ny[sel[i]])] = y[,sel[i]]
	
	if(length(nn)>0) {
		yinf[,((length(nx)+1):length(nm))] = y[,!(ny%in%nx)] #Remaining columns in y
	}
	row.names(yinf) = row.names(y)
	names(yinf) = nm
	
	return (list(x=xinf, y =yinf))
}