#' Heterogeneity-constrained random resampling (HCR) 
#' 
#' Returns a set of indices of the original data set that maximizes the mean and minimizes the variance of the distances between pairs of plot records.
#' 
#' @param d An object of class \code{\link{dist}} containing the distance values between pairs of sites (plot records).
#' @param nout The number of sites (plot records) to be chosen among those available in \code{d}.
#' @param nsampl The number of resampling trials to be compared.
#'
#' @details
#' Many subsets of the input data are selected randomly. These subsets are sorted by decreasing mean dissimilarity between pairs of plot records, and then sorted again by increasing variance of these dissimilarities. Ranks from both sortings are summed for each subset, and the subset with the lowest summed rank is considered as the most representative.
#' 
#' @returns
#' Returns a vector containing the indices of the selected sites (plot records) to be used for sub-setting the original table.
#' 
#' @references 
#' Lengyel, A., Chytry, M., Tichy, L. (2011) Heterogeneity-constrained random resampling of phytosociological databases. Journal of Vegetation Science 22: 175-183.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' 
#' @seealso \code{\link[vegan]{decostand}}
#' 
#' @export
#'
#' @examples
#' ## Loads data (38 columns and 33 species)
#' data(wetland)
#' dim(wetland)
#' 
#' ## Constructs the chord distance matrix
#' wetland.chord <-dist(as.data.frame(sweep(as.matrix(wetland), 1, 
#'                                          sqrt(rowSums(as.matrix(wetland)^2)), "/")))
#' 
#' ## Performs HCR resampling. Returns indices of objects
#' sel <- hcr(wetland.chord, nout=20, nsampl=1000)
#' 
#' ## Prints the names of the plot records
#' print(row.names(wetland)[sel])
#' 
#' ## Subset the original distance matrix
#' sel.chord <- as.dist(as.matrix(wetland.chord)[sel,sel])
hcr<-function(d, nout, nsampl=1000){
	dmat = as.matrix(d)
	nplots = nrow(dmat)
	if(nout > nplots) stop("nout cannot be larger than the number of plots")
	meand = numeric(nsampl)
	vard = numeric(nsampl)
	sel<-matrix(FALSE,nrow=nout, ncol=nsampl)
	for(s in 1:nsampl) {
		sel[,s] <- sample(x=nplots,size=nout)
		dvec = as.vector(as.dist(dmat[sel[,s],sel[,s]]))
		meand[s] = mean(dvec)
		vard[s] = var(dvec)
	}
	#print(head(sort(meand, decreasing=TRUE)))
	rankdecmean = rank(-meand)
	#print(head(sort(vard)))	
	rankincvar = rank(vard)
	rankfinal = rank(rankdecmean+rankincvar)
	#print(cbind(meand, rankdecmean, vard, rankincvar,rankfinal))
	finalsel = sort(sel[,which.min(rankfinal)])
	return(finalsel)
}