#' Classifies vegetation communities
#' 
#' Classifies vegetation communities into a previous fuzzy or hard classification.
#'
#' @param y An object of class \code{\link{vegclust}} that represents a previous knowledge.
#' @param x Community data to be classified, in form of a site by species matrix (if the vegclust object is in \code{raw} mode) or a data frame containing the distances between the new sites in rows and the old sites in columns (if the \code{\link{vegclust}} object is in \code{distance} mode).
#'
#' @details
#' This function uses the classification model specified in \code{y} to classify the communities (rows) in \code{x}. When vegclust is in \code{raw} mode, the function calls first to \code{\link{conformveg}} in order to cope with different sets of species. See the help of \code{\link{as.vegclust}} to see an example of \code{vegclass} with distance matrices.
#' 
#' @returns Returns an object of type \code{vegclass} with the following items:
#' \itemize{
#' \item{\code{method}: The clustering model used in \code{y}}
#' \item{\code{m}: The fuzziness exponent in \code{y}}
#' \item{\code{dnoise}:The distance to the noise cluster used for noise clustering (models NC, NCdd, HNC, HNCdd). This is set to \code{NULL} for other models.}
#' \item{\code{eta}: The reference distance vector used for possibilistic clustering (models PCM and PCMdd). This is set to \code{NULL} for other models.}
#' \item{\code{memb}: The fuzzy membership matrix.}
#' \item{\code{dist2clusters}: The matrix of object distances to cluster centers.}
#' }
#' 
#' @references
#' \enc{Davé}{Dave}, R. N. and R. Krishnapuram (1997) Robust clustering methods: a unified view. IEEE Transactions on Fuzzy Systems 5, 270-293.
#' 
#' Bezdek, J. C. (1981) Pattern recognition with fuzzy objective functions. Plenum Press, New York.
#' 
#' Krishnapuram, R. and J. M. Keller. (1993) A possibilistic approach to clustering. IEEE transactions on fuzzy systems 1, 98-110.
#' 
#' De \enc{Cáceres}{Caceres}, M., Font, X, Oliva, F. (2010) The management of numerical vegetation classifications with fuzzy clustering methods [Related software]. Journal of Vegetation Science 21 (6): 1138-1151.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF.
#' 
#' @seealso \code{\link{vegclust}}, \code{\link{as.vegclust}}, \code{\link{kmeans}}, \code{\link{conformveg}}
#' 
#' @export
#'
#' @examples
#' ## Loads data (38 columns and 33 species)
#' data(wetland)
#' 
#' ## This equals the chord transformation 
#' wetland.chord <- as.data.frame(sweep(as.matrix(wetland), 1, 
#'                                      sqrt(rowSums(as.matrix(wetland)^2)), "/"))
#' 
#' ## Splits wetland data into two matrices of 30x27 and 11x22
#' wetland.30 <- wetland.chord[1:30,]
#' wetland.30 <- wetland.30[,colSums(wetland.30)>0]
#' dim(wetland.30)
#' wetland.11 <- wetland.chord[31:41,]
#' wetland.11 <- wetland.11[,colSums(wetland.11)>0] 
#' dim(wetland.11)
#' 
#' ## Create noise clustering with 3 clusters from the data set with 30 sites. 
#' wetland.30.nc <- vegclust(wetland.30, mobileCenters=3, m = 1.2, dnoise=0.75,
#'                           method="NC", nstart=10)
#' 
#' ## Cardinality of fuzzy clusters (i.e., the number of objects belonging to)
#' wetland.30.nc$size
#' 
#' ## Classifies the second set of sites according to the clustering of the first set
#' wetland.11.nc <- vegclass(wetland.30.nc, wetland.11)
#' 
#' ## Fuzzy membership matrix
#' wetland.11.nc$memb
#' 
#' ## Obtains hard membership vector, with 'N' for objects that are unclassified
#' defuzzify(wetland.11.nc$memb)$cluster
vegclass<-function(y, x) {
	if(!inherits(y,"vegclust")) stop("y must be a vegclust object")
	
	if(y$mode=="raw"){
	  if(!is.null(y$fixedCenters)) centers = rbind(y$mobileCenters,y$fixedCenters)
	  else centers = y$mobileCenters
	  if(length(colnames(x))!=length(colnames(centers)) || sum(colnames(x)==colnames(centers))<ncol(x)) {
	    c = conformveg(x,centers)
	    x = as.matrix(c$x)
	    centers = as.matrix(c$y)
	  } else {
	    x = as.matrix(x)
	    centers = as.matrix(centers)
	  }
	  k = nrow(centers)
	} else { 
	  memb = as.matrix(y$memb)
	  x = as.matrix(x) #x contains the distance from new objects to old ones
	  d2cl = as.matrix(y$dist2clusters)
	  k = ncol(d2cl)
	}
	
	m = y$m
  dnoise = y$dnoise
  eta = y$eta
  method = y$method
	n = nrow(x)
		

  if(method=="NC"||method=="NCdd"||method=="HNC"||method=="HNCdd") {
   	u = matrix(0,nrow=n,ncol=(k+1))
  } else {
   	u = matrix(0,nrow=n,ncol=k)
  }
	

	dist2cent = matrix(0,nrow=nrow(x),ncol=k)

	#1. compute distance to centers (fixed and mobile)
	if(y$mode=="raw") {
	  	for(i in 1:k) {
		   dist2cent[,i] = sqrt(rowSums(sweep(x,2,centers[i,],"-")^2))
		}
	} else { #distance mode
	  if(method=="KM"||method=="PCM"||method=="NC"|| method=="HNC"||method=="FCM") {     
  		for(i in 1:k) {
	  		vargeom = sum((memb[,i]^m) %*% (d2cl[,i]^2))/sum(memb[,i]^m)
		  	for(j in 1:n) {
	  			dist2cent[j,i] = sqrt((sum((memb[,i]^m)*(x[j,]^2))/sum(memb[,i]^m))-vargeom)
			  }
		  }
	  } else { # For medoid methods, medoids are stored in mobileCenters and fixedCenters 
       med = c(y$mobileCenters,y$fixedCenters)
	     for(i in 1:k) {
	       dist2cent[,i] = x[,med[i]] 
	      }
	  }
	}
	
	  #2. compute membership to centroids
	if (y$method=="KM"||y$method=="KMdd") {
	  minC<-apply(dist2cent,1,which.min)
	  u[,] = 0
	  for(k in 1:length(minC)) u[k,minC[k]] = 1.0
	} else if(y$method=="NC") {
	  d2cm2<-cbind(dist2cent,dnoise)^2
	  for(k in 1:ncol(d2cm2)) {
	    a<-sweep(d2cm2,1,d2cm2[,k],"/")
	    u[,k] = 1/rowSums(a^(-1/(m-1)))
	  }
	  u[d2cm2==0]=1
	} else if(y$method=="HNC"||y$method=="HNCdd") {
	  d2cm<-cbind(dist2cent,dnoise)
	  u[,] = 0
	  minC<-apply(d2cm,1,which.min)
	  for(k in 1:length(minC)) {
	    u[k,minC[k]] = 1.0
	  }
	} else if(y$method=="NCdd") {
	  d2cm<-cbind(dist2cent,dnoise)
	  for(k in 1:ncol(d2cm)) {
	    a<-sweep(d2cm,1,d2cm[,k],"/")
	    u[,k] = 1/rowSums(a^(-1/(m-1)))
	  }
	  u[d2cm==0]=1
	} else if (y$method=="FCM") {
	  d2cm2<-dist2cent^2
	  for(k in 1:ncol(dist2cent)) {
	    a<-sweep(d2cm2,1,d2cm2[,k],"/")
	    u[,k] = 1/rowSums(a^(-1/(m-1)))
	  }
	  u[dist2cent ==0]=1
	} else if (y$method=="FCMdd") {
	  d2cm<-dist2cent
	  for(k in 1:ncol(d2cm)) {
	    a<-sweep(d2cm,1,d2cm[,k],"/")
	    u[,k] = 1/rowSums(a^(-1/(m-1)))
	  }
	  u[dist2cent ==0]=1
	} else if (y$method=="PCM") {
	  for(k in 1:ncol(dist2cent)) u[,k] = 1/(1+((dist2cent[,k]^2)/eta[k])^(1/(m-1)))
	  u[dist2cent==0]=1
	} else if (y$method=="PCMdd") {
	  for(k in 1:ncol(dist2cent)) u[,k] = 1/(1+((dist2cent[,k])/eta[k])^(1/(m-1)))
	  u[dist2cent==0]=1
	} 	 	
	
   #Prepare output
   u = as.data.frame(u)   
   dist2cent = as.data.frame(dist2cent)   
   if(ncol(u)==ncol(y$memb)) names(u) = names(y$memb)
   else {
   	names(u)[1:ncol(y$memb)] = names(y$memb)
   	names(u)[ncol(y$memb)+1] = "N"
   }
   names(dist2cent) = names(y$dist2clusters)   
	rownames(u) = rownames(x)
	rownames(dist2cent) = rownames(x)
	
   res = list(method = y$method, m =y$m, dnoise = y$dnoise, eta= y$eta, memb=u,dist2clusters=dist2cent)
   class(res)<-"vegclass"
	return(res)
		
}