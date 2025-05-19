#' Relates two clustering level results
#' 
#' Analyzes how lower level clusters are assigned into upper level ones. The analysis is made for several number of clusters.
#' 
#' @param lower A list of objects of type \code{\link{vegclust}} or \code{\link{vegclass}} that represent classifications at a finer level of resolution.
#' @param upper A list of objects of type \code{\link{vegclust}} or \code{\link{vegclass}} that represent classifications at a broader level of resolution.
#' @param defuzzify A logical flag used to indicate whether the result of calling \code{\link{crossmemb}} should be deffuzified.
#' @param excludeFixed A logical used to indicate whether fixed clusters should be excluded from the comparison of levels.
#' @param verbose A flag used to ask for extra screen output.
#' @param ... Additional parameters for function \code{\link{defuzzify}}.
#'
#' @details
#' For each pair of \code{vegclust} (or \code{vegclass}) objects in \code{upper} and \code{lower}, the function calls function \code{\link{crossmemb}} and then, if asked, deffuzifies the resulting memberships (by calling function \code{\link{defuzzify}}) and several quantities are calculated (see 'value' section).
#' 
#' @returns  A list with several data frames (see below). In each of them, the rows are items of \code{upper} and columns are items of \code{lower}. The names of rows and columns are the number of clusters of each \code{\link{vegclust}} (or \code{vegclass}) object.
#' \itemize{
#' \item{\code{nnoise}: The number of low level clusters that are assigned to the Noise class (for \code{upper} objects using Noise clustering). }
#' \item{\code{maxnoise}: The maximum membership value of low level clusters to the Noise class (for \code{upper} objects using Noise clustering). }
#' \item{\code{minmaxall}: The minimum value (across upper level clusters) of the maximum membership value observed among the lower level clusters. }
#' \item{\code{minallsize}: The minimum value (across upper level clusters) of the sum of membership values across lower level clusters. }
#' \item{\code{empty}: The number of upper level clusters (mobile or fixed) that do not have any member among the lower level clusters. }
#' }
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' @seealso \code{\link{vegclust}}, \code{\link{vegclass}}, \code{\link{defuzzify}}
#' 
#' @export
#'
#' @examples
#' ## Loads data  
#' data(wetland)
#' 
#' ## This equals the chord transformation 
#' wetland.chord <- as.data.frame(sweep(as.matrix(wetland), 1,
#'                                      sqrt(rowSums(as.matrix(wetland)^2)), "/"))
#' 
#' ## Create noise clustering from hierarchical clustering at different number of cluster
#' wetland.hc <- hclust(dist(wetland.chord),method="ward") 
#' wetland.nc1 <- hier.vegclust(wetland.chord, wetland.hc, cmin=2, cmax=6, m = 1.2, 
#'                              dnoise=0.75, method="NC")
#' wetland.nc2 <- hier.vegclust(wetland.chord, wetland.hc, cmin=2, cmax=4, m = 1.2, 
#'                              dnoise=0.85, method="NC")
#' 
#' ## Studies the assignment of levels
#' relate.levels(wetland.nc1, wetland.nc2, method="cut")
#' 
relate.levels<-function(lower, upper, defuzzify=FALSE, excludeFixed = FALSE, verbose=FALSE, ...) {
  minupperclasses =999999
  maxupperclasses =0
	minlowerclasses=999999
  	maxlowerclasses =0
	for(i in 1:length(upper)) {
		if(!is.null(upper[[i]])) {
			nclasses = length(names(upper[[i]]$dist2clusters))
			maxupperclasses = max(maxupperclasses,nclasses)
			minupperclasses = min(minupperclasses,nclasses)
		}
	}
	for(i in 1:length(lower)) {
		if(!is.null(lower[[i]])) {
			nclasses = length(names(lower[[i]]$dist2clusters))
			minlowerclasses = min(minlowerclasses,nclasses)
			maxlowerclasses = max(maxlowerclasses,nclasses)
		}
	}
	#Prepare output matrices
	maxnoise = data.frame(matrix(NA,length(minupperclasses:maxupperclasses), length(minlowerclasses:maxlowerclasses)))
	row.names(maxnoise) = minupperclasses:maxupperclasses
	names(maxnoise) = minlowerclasses:maxlowerclasses
	minmaxall = maxnoise
	nnoise = maxnoise
	minallsize = maxnoise
	empty = maxnoise
	
	#Loop over lower clustering
	for(j in 1:length(lower)) {
		if(!is.null(lower[[j]])) {
			nlowerclasses = length(names(lower[[j]]$dist2clusters))
			colind<-nlowerclasses-minlowerclasses+1
			if(excludeFixed) {
			  nfixed = sum(substr(names(lower[[j]]$dist2clusters),1,1)=="F")
        nlowerclasses = nlowerclasses - nfixed
        colind<-nlowerclasses-minlowerclasses+1 + nfixed
			}
			if(verbose) cat(paste("Number of lower classes:", nlowerclasses,"\n"))
			#Loop over upper clustering
			for(i in 1:length(upper)) {
				if(!is.null(upper[[i]])){
					if(verbose) cat(".")
					nupperclasses = length(names(upper[[i]]$dist2clusters))
					rowind<-nupperclasses-minupperclasses+1
					if(excludeFixed) {
            nfixed = sum(substr(names(upper[[i]]$dist2clusters),1,1)=="F")
            nupperclasses = nupperclasses - nfixed
            rowind<-nupperclasses-minupperclasses+1 + nfixed
					}
					memb<-crossmemb(lower[[j]], upper[[i]])
					if(lower[[j]]$method=="NC") memb = memb[-nrow(memb),]
					if(defuzzify) memb <- defuzzify(memb, ...)$memb
					minmaxall[rowind,colind] = min(apply(memb[1:nlowerclasses,1:nupperclasses],2,max))
					minallsize[rowind,colind] = min(apply(memb[1:nlowerclasses,1:nupperclasses],2,sum))
					empty[rowind,colind] = sum(colSums(defuzzify(memb[1:nlowerclasses,1:nupperclasses], method="max")$memb)==0)
					#print(empty[rowind,colind])
					if(upper[[i]]$method=="NC") {
						maxnoise[rowind,colind]=max(memb[1:nlowerclasses,ncol(memb)])
						nplots = colSums(defuzzify(lower[[j]], method="cut", alpha=0.5)$memb)
						nnoise[rowind,colind]=sum(memb[1:nlowerclasses,ncol(memb)]>0.5)
					 	#print(nnoise[rowind,colind])
					}
				}
			}
		}
		if(verbose) cat("done.\n")
	}	

	return(list(nnoise =nnoise, maxnoise=maxnoise, minallsize=minallsize, minmaxall = minmaxall, empty= empty))
}