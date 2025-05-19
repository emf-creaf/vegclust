#' Vegetation clustering methods
#' 
#' Performs hard or fuzzy clustering of vegetation data
#'
#' @param x Community data. A site-by-species matrix or data frame (for \code{vegclust}) or a site-by-site dissimilarity matrix or \code{\link{dist}} object (for \code{vegclustdist}).
#' @param mobileCenters A number, a vector of seeds, or coordinates for mobile clusters.
#' @param fixedCenters A matrix or data frame with coordinates for fixed (non-mobile) clusters.
#' @param method A clustering model. Current accepted models are: 
#' \itemize{
#'   \item{\code{"KM"}: K-means or hard c-means (MacQueen 1967)}
#'   \item{\code{"KMdd"}: Hard c-medoids (Krishnapuram et al. 1999)}
#'   \item{\code{"FCM"}: Fuzzy c-means (Bezdek 1981)}
#'   \item{\code{"FCMdd"}: Fuzzy c-medoids (Krishnapuram et al. 1999)}
#'   \item{\code{"NC"}: Noise clustering (Dave and Krishnapuram 1997)}
#'   \item{\code{"NCdd"}: Noise clustering with medoids}
#'   \item{\code{"HNC"}: Hard noise clustering}
#'   \item{\code{"HNCdd"}: Hard noise clustering with medoids}
#'   \item{\code{"PCM"}: Possibilistic c-means (Krishnapuram and Keller 1993)}
#'   \item{\code{"PCMdd"}: Possibilistic c-medoids}
#' }
#' @param m The fuzziness exponent to be used (this is relevant for all models except for kmeans)
#' @param dnoise The distance to the noise cluster, relevant for noise clustering (NC). 
#' @param eta A vector of reference distances, relevant for possibilistic C-means (PCM). 
#' @param alpha Threshold used to stop iterations. The maximum difference in the membership matrix of the current vs. the previous iteration will be compared to this value.
#' @param iter.max The maximum number of iterations allowed.
#' @param nstart If \code{mobileCenters} or \code{mobileMemb} is a number, how many random sets should be chosen?
#' @param maxminJ When random starts are used, these will stop if at least \code{maxminJ} runs ended up in the same functional value.
#' @param seeds If \code{mobileCenters} or \code{mobileMemb} is a number, a vector indicating which objects are potential initial centers. If \code{NULL} all objects are valid seeds.
#' @param verbose Flag to print extra output.
#' 
#' @details 
#' Functions \code{vegclust} and \code{vegclustdist} try to generalize the \code{\link{kmeans}} function in \code{stats} in three ways. 
#' 
#' Firstly, they allows different clustering models. Clustering models can be divided in (a) fuzzy or hard; (b) centroid-based or medoid-based; (c) Partitioning (KM and FCM family), noise clustering (NC family), and possibilistic clustering (PCM and PCMdd). The reader should refer to the original publications to better understand the differences between models. 
#' 
#' Secondly, users can specify fixed clusters (that is, centroids that do not change their positions during iterations). Fixed clusters are intended to be used when some clusters were previously defined and new data has been collected. One may allow some of these new data points to form new clusters, while some other points will be assigned to the original clusters. In the case of models with cluster repulsion (such as KM, FCM or NC) the new (mobile) clusters are not allowed to 'push' the fixed ones. As a result, mobile clusters will occupy new regions of the reference space. 
#' 
#' Thirdly, \code{vegclustdist} implements the distance-based equivalent of \code{vegclust}. The results of \code{vegclust} and \code{vegclustdist} will be the same (if seeds are equal) if the distance matrix is calculated using the Euclidean distance (see function \code{\link{dist}}). Otherwise, the equivalence holds by resorting on principal coordinates analysis.
#' 
#' Note that all data frames or matrices used as input of \code{vegclust} should be defined on the same space of species (see \code{\link{conformveg}}). Unlike \code{\link{kmeans}}, which allows different specific algorithms, here updates of prototypes (centroids or medoids) are done after all objects have been reassigned (Forgy 1965). In order to obtain hard cluster definitions, users can apply the function \code{\link{defuzzify}} to the \code{vegclust} object.
#' 
#' @returns Returns an object of type \code{vegclust} with the following items:
#' \itemize{
#' \item{\code{mode}: \code{raw} for function \code{vegclust} and \code{dist} for function \code{vegclustdist}.}
#' \item{\code{method}: The clustering model used}
#' \item{\code{m}: The fuzziness exponent used (\code{m=1} in case of kmeans)}
#' \item{\code{dnoise}: The distance to the noise cluster used for noise clustering (NC, HNC, NCdd or HNCdd). This is set to \code{NULL} for other models.}
#' \item{\code{eta}: The reference distance vector used for possibilistic clustering (PCM or PCMdd). This is set to \code{NULL} for other models.}
#' \item{\code{memb}: The fuzzy membership matrix. Columns starting with "M" indicate mobile clusters, whereas columns starting with "F" indicate fixed clusters.}
#' \item{\code{mobileCenters}: If \code{vegclust} is used, this contains a data frame with the coordinates of the mobile centers (centroids or medoids). If \code{vegclustdist} is used, it will contain the indices of mobile medoids for models KMdd, FCMdd, HNCdd, NCdd and PCMdd; or \code{NULL} otherwise.}
#' \item{\code{fixedCenters}: If \code{vegclust} is used, this contains a data frame with the coordinates of the fixed centers (centroids or medoids). If \code{vegclustdist} is used, it will contain the indices of fixed medoids for models KMdd, FCMdd, HNCdd, NCdd and PCMdd; or \code{NULL} otherwise.}
#' \item{\code{dist2clusters}: The matrix of object distances to cluster centers. Columns starting with "M" indicate mobile clusters, whereas columns starting with "F" indicate fixed clusters.}
#' \item{\code{withinss}: In the case of methods KM, FCM, NC, PCM and HNC it contains the within-cluster sum of squares for each cluster (squared distances to cluster center weighted by membership). In the case of methods KMdd, FCMdd, NCdd, HNCdd and PCMdd it contains the sum of distances to each cluster (weighted by membership).}
#' \item{\code{size}: The number of objects belonging to each cluster. In case of fuzzy clusters the sum of memberships is given.}
#' \item{\code{functional}: The objective function value (the minimum value attained after all iterations).}
#' }
#' 
#' @references
#' Forgy, E. W. (1965) Cluster analysis of multivariate data: efficiency vs interpretability of classifications. Biometrics 21, 768-769.
#' 
#' MacQueen, J. (1967) Some methods for classification and analysis of multivariate observations. In Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability, eds L. M. Le Cam and J. Neyman, 1, pp. 281-297. Berkeley, CA: University of California Press.
#' 
#' \enc{Davé}{Dave}, R. N. and R. Krishnapuram (1997) Robust clustering methods: a unified view. IEEE Transactions on Fuzzy Systems 5, 270-293.
#' 
#' Bezdek, J. C. (1981) Pattern recognition with fuzzy objective functions. Plenum Press, New York.
#' 
#' Krishnapuram, R., Joshi, A., & Yi, L. (1999). A Fuzzy relative of the k-medoids algorithm with application to web document and snippet clustering. IEEE International Fuzzy Systems (pp. 1281–1286). 
#' 
#' Krishnapuram, R. and J. M. Keller. (1993) A possibilistic approach to clustering. IEEE transactions on fuzzy systems 1, 98-110.
#' 
#' De \enc{Cáceres}{Caceres}, M., Font, X, Oliva, F. (2010) The management of numerical vegetation classifications with fuzzy clustering methods. Journal of Vegetation Science 21 (6): 1138-1151.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' 
#' @seealso \code{\link{hier.vegclust}},\code{\link{incr.vegclust}},\code{\link{kmeans}},\code{\link{vegclass}},\code{\link{defuzzify}},\code{\link{clustvar}}
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
#' ## Create noise clustering with 3 clusters. Perform 10 starts from random seeds 
#' ## and keep the best solution
#' wetland.nc <- vegclust(wetland.chord, mobileCenters=3, m = 1.2, dnoise=0.75, 
#'                        method="NC", nstart=10)
#' 
#' ## Fuzzy membership matrix
#' wetland.nc$memb
#' 
#' ## Cardinality of fuzzy clusters (i.e., the number of objects belonging to each cluster)
#' wetland.nc$size
#' 
#' ## Obtains hard membership vector, with 'N' for objects that are unclassified
#' defuzzify(wetland.nc$memb)$cluster
#' 
#' ## The same result is obtained with a matrix of chord distances
#' wetland.d <- dist(wetland.chord)
#' wetland.d.nc <- vegclustdist(wetland.d, mobileMemb=3, m = 1.2, dnoise=0.75, 
#'                              method="NC", nstart=10)
vegclust <-
function(x,mobileCenters, fixedCenters=NULL, method="NC", m=2,dnoise=NULL, eta = NULL, alpha=0.001, iter.max=100, nstart=1, maxminJ=10, seeds = NULL, verbose=FALSE) {

#One run of vegclust   
vegclustone <-
function(x,mobileCenters, fixedCenters=NULL, method="NC", m=2,dnoise=NULL, eta = NULL, alpha=0.001, iter.max=100) {
   METHODS <- c("KM", "FCM", "PCM","NC","HNC" ,"KMdd","NCdd", "HNCdd", "FCMdd", "PCMdd")
   method <- match.arg(method, METHODS)
   if(method=="KM"||method=="KMdd") {
   	m=1.0
   	dnoise=NULL
   	eta=NULL
   }
   else if(method=="FCM"||method=="FCMdd") {
   	dnoise=NULL
   	eta=NULL
   }
   else if(method=="NC"||method=="NCdd") {
     if(is.null(dnoise)) stop("Must provide a value for dnoise")
     eta = NULL
   }
   else if(method=="HNC"||method=="HNCdd") {
     if(is.null(dnoise)) stop("Must provide a value for dnoise")
     eta = NULL
     m=1.0
   }
   else if(method=="PCM"||method=="PCMdd") {
   	if(is.null(eta)) stop("Must provide a vector of values for eta")
   	dnoise = NULL
   }
   
	x = as.matrix(x)
  #If medoid methods are used, calculate the (Euclidean) distance matrix
	if(method=="KMdd"||method=="FCMdd"||method=="NCdd"||method=="PCMdd"||method=="HNCdd") {
		d = as.matrix(dist(x))
	}
  #Number of objects
	n = nrow(x)
	
	#Sets the starting mobile centers (can be centroids or medoids)
	if(is.data.frame(mobileCenters) || is.matrix(mobileCenters)) {
		mobileCenters = as.matrix(mobileCenters)
		if(ncol(mobileCenters)!=ncol(x)) {
			stop("The number and identity of species for mobile centers must be the same as the number and identity of species for x")
		}		
	}
	else if(is.vector(mobileCenters) && length(mobileCenters)==1 && is.numeric(mobileCenters)) {
		if(mobileCenters==1) mobileCenters = t(as.matrix(x[sample(n,mobileCenters),]))
		else mobileCenters = as.matrix(x[sample(n,mobileCenters),])
	}
	else if(is.vector(mobileCenters) && is.numeric(mobileCenters)) {
		mobileCenters = as.matrix(x[mobileCenters,])
	} 
	if(!is.matrix(mobileCenters)) {
		stop("Provide a number, a vector of seeds, or coordinates for mobile centers")
	}	
	kMov = nrow(mobileCenters)
	rownames(mobileCenters)<-c(1:kMov)
	
	#Sets the fixed centers (can be centroids or medoids)
	if(!is.null(fixedCenters)) {
		if(is.data.frame(fixedCenters)) {
			fixedCenters = as.matrix(fixedCenters)
			kFix = nrow(fixedCenters)
		}
		else if(!is.matrix(fixedCenters)) {
			stop("Fixed centers must be specified as a matrix or a data frame")
		}	
		else {kFix = nrow(fixedCenters)}
		if(ncol(fixedCenters)!=ncol(x)) {
			stop("The number and identity of species for fixed centers must be the same as the number and identity of species for x")
		}
	} else {
		kFix = 0
	}

   if((method=="PCM"||method=="PCMdd") && length(eta)!=(kMov+kFix)) stop("Vector of reference distances (eta) must have a length equal to the number of clusters")
  	dist2cent = matrix(0,nrow=n,ncol=(kMov+kFix))

   #Add an extra (noise) column for NC-related methods
   if(method=="NC"||method=="NCdd"||method=="HNC"||method=="HNCdd") {
   		u = matrix(0,nrow=n,ncol=(kMov+kFix+1))
  		uPrev = matrix(0,nrow=n,ncol=(kMov+kFix+1))
  	} else {
   		u = matrix(0,nrow=n,ncol=(kMov+kFix))
  		uPrev = matrix(0,nrow=n,ncol=(kMov+kFix))
  	}

	  #1. compute distance to mobile and fixed centers 
   	for(k in 1:kMov) {
   		dist2cent[,k] = sqrt(rowSums(sweep(x,2,mobileCenters[k,],"-")^2))
   	}
  	if(kFix==1) {
  		dist2cent[,1+kMov] = sqrt(rowSums(sweep(x,2, as.vector(fixedCenters),"-")^2))
  	}  			
  	else if(kFix>1) {
  		for(k in 1:kFix) {
  			dist2cent[,k+kMov] = sqrt(rowSums(sweep(x,2,fixedCenters[k,],"-")^2))
  		}
  	}  			

	continue = TRUE
	iter = 1
   #iterates until no change in memberships
   while(continue) {
	  #1. compute membership to centers (centroids or medoids)
   	if (method=="KM"||method=="KMdd") {
   	    minC<-apply(dist2cent,1,which.min)
        u[,] = 0
   	    for(k in 1:length(minC)) u[k,minC[k]] = 1.0
		} else if(method=="NC") {
   			d2cm2<-cbind(dist2cent,dnoise)^2
   			for(k in 1:ncol(d2cm2)) {
   				a<-sweep(d2cm2,1,d2cm2[,k],"/")
   				u[,k] = 1/rowSums(a^(-1/(m-1)))
   			}
			u[d2cm2==0]=1
		} else if(method=="HNC"||method=="HNCdd") {
		    d2cm<-cbind(dist2cent,dnoise)
		    u[,] = 0
		    minC<-apply(d2cm,1,which.min)
		    for(k in 1:length(minC)) {
          u[k,minC[k]] = 1.0
		    }
		} else if(method=="NCdd") {
   			d2cm<-cbind(dist2cent,dnoise)
   			for(k in 1:ncol(d2cm)) {
   				a<-sweep(d2cm,1,d2cm[,k],"/")
   				u[,k] = 1/rowSums(a^(-1/(m-1)))
   			}
			u[d2cm==0]=1
		} else if (method=="FCM") {
   			d2cm2<-dist2cent^2
   			for(k in 1:ncol(dist2cent)) {
   				a<-sweep(d2cm2,1,d2cm2[,k],"/")
   				u[,k] = 1/rowSums(a^(-1/(m-1)))
   			}
			u[dist2cent ==0]=1
		} else if (method=="FCMdd") {
   			d2cm<-dist2cent
   			for(k in 1:ncol(d2cm)) {
   				a<-sweep(d2cm,1,d2cm[,k],"/")
   				u[,k] = 1/rowSums(a^(-1/(m-1)))
   			}
			u[dist2cent ==0]=1
		} else if (method=="PCM") {
			for(k in 1:ncol(dist2cent)) u[,k] = 1/(1+((dist2cent[,k]^2)/eta[k])^(1/(m-1)))
			u[dist2cent==0]=1
		} else if (method=="PCMdd") {
			for(k in 1:ncol(dist2cent)) u[,k] = 1/(1+((dist2cent[,k])/eta[k])^(1/(m-1)))
			u[dist2cent==0]=1
		} 	 	
   		#Check for stopping
   		if(iter>2) {
   			continue = (max(abs(u-uPrev))>alpha) && (iter<=iter.max) && (max(abs(u-uPrev2))>alpha)
   		}
   	
   	if(continue) {
	   	#2. update mobile centers (centroids or medoids) and distances
	   	if(method=="KM"||method=="PCM"||method=="NC"|| method=="HNC"||method=="FCM") {
		   	for(k in 1:kMov) {
		   		mobileCenters[k,]=((u[,k]^m)%*%x)/sum(u[,k]^m)	   		
		   		dist2cent[,k] = sqrt(rowSums(sweep(x,2,mobileCenters[k,],"-")^2))
   			}
	   	} else { 
	   		for(k in 1:kMov) {
	   			#Determine medoid
	   			med = which.min((u[,k]^m)%*%d)
	   			#Update mobile medoids and distances
	   			mobileCenters[k,] = x[med,]
	   			dist2cent[,k] = d[,med]
	   		}
	   	}
	   	uPrev2 = uPrev
	   	uPrev = u	   	
	   	iter=iter+1
	   	if(verbose) cat(".")
   	}
   }
  	if(method=="FCM" ||  method=="KM") functional = sum((dist2cent^2)*(u^m))
  	else if(method=="NC"||method=="HNC") functional = sum((dist2cent^2)*(u[,-(kMov+kFix+1)]^m))+sum(dnoise^2*u[,kMov+kFix+1]^m)
  	else if(method=="FCMdd"||  method=="KMdd") functional = sum(dist2cent*(u^m))
  	else if(method=="NCdd"||method=="HNCdd") functional = sum(dist2cent*(u[,-(kMov+kFix+1)]^m))+sum(dnoise*u[,kMov+kFix+1]^m)
  	else if(method=="PCM") {
  		functional = 0
  		for(k in 1:(kMov+kFix)) functional = functional+sum((dist2cent[,k]^2)*(u[,k]^m))+sum(eta[k]*(1-u[,k])^m)
  	} else if(method=="PCMdd") {
  		functional = 0
  		for(k in 1:(kMov+kFix)) functional = functional+sum(dist2cent[,k]*(u[,k]^m))+sum(eta[k]*(1-u[,k])^m)
  	} 
   if(verbose) cat(paste("\nIterations:", iter,"Functional: ", functional,"\n"))
   #Prepare output
   mobileCenters = as.data.frame(mobileCenters)
   u = as.data.frame(u)   
   dist2cent = as.data.frame(dist2cent)   
	for(k in 1:kMov) {
		rownames(mobileCenters)[k] = paste("M",k,sep="")
		colnames(u)[k] = paste("M",k,sep="")
		colnames(dist2cent)[k] = paste("M",k,sep="")
	}
	if(kFix>1) {
		for(k in (kMov+1):(kMov+kFix)) {
			colnames(u)[k] = paste("F",k,sep="")
			colnames(dist2cent)[k] = paste("F",k,sep="")
		}
	}
	if(method=="NC"||method=="NCdd"||method=="HNC"||method=="HNCdd") names(u)[kMov+kFix+1] = "N"
	rownames(u) = rownames(x)
	rownames(dist2cent) = rownames(x)
	size = colSums(u[,1:(kMov+kFix), drop=FALSE])
   if(method=="NC"||method=="FCM"||method=="KM"||method=="PCM"||method=="HNC") withinss = colSums((dist2cent^2)*(u[,1:(kMov+kFix), drop=FALSE]^m))
   else withinss = colSums((dist2cent)*(u[,1:(kMov+kFix), drop=FALSE]^m))
   res = list(mode="raw", method=method, m = m, dnoise = dnoise,eta = eta, memb=u,mobileCenters=mobileCenters, fixedCenters=fixedCenters, dist2clusters=dist2cent, withinss = withinss, size=size, functional=functional, iter=iter)
   class(res)<-"vegclust"
	return(res)
}

	if(is.null(seeds)) seeds = 1:nrow(x)
	#print(seeds)
   #If mobileCenters is a number and nstart>1 perform different random starts
	if(is.vector(mobileCenters) && length(mobileCenters)==1 && is.numeric(mobileCenters)) {
	   bestRun = vegclustone(x,mobileCenters=x[sample(seeds,mobileCenters),], fixedCenters, method, m,dnoise, eta, alpha, iter.max)
		if(nstart>1) {
			minJ = 0
			i = 2
			while(i<=nstart) {
				run = vegclustone(x,mobileCenters=x[sample(seeds,mobileCenters),], fixedCenters, method, m,dnoise,eta, alpha, iter.max)
				if(run$functional<=bestRun$functional) {
					bestRun = run
					if(max(abs(run$memb-bestRun$memb))<alpha) {
						minJ = minJ+1
					} else {
						minJ = 1
					}
				}
				if(minJ==maxminJ) {
					i=nstart
					if(verbose) print("Maximum number of minimum J reached. Stopping.")
				}
				i=i+1
			}
		}		
		return(bestRun)
	} else { #Perform a single run
		return(vegclustone(x,mobileCenters, fixedCenters, method, m,dnoise, eta, alpha, iter.max))
	}
}


#' @rdname vegclust
#' @export
#' @param mobileMemb A number, a vector of seeds, or starting memberships for mobile clusters.
#' @param fixedDistToCenters A matrix or data frame with the distances to fixed cluster centers.
vegclustdist <-
  function(x,mobileMemb, fixedDistToCenters=NULL, method="NC", m=2,dnoise=NULL, eta = NULL, alpha=0.001, iter.max=100, nstart=1, seeds = NULL, verbose=FALSE) {
    
    #One run of vegclustdist   
    vegclustonedist <-
      function(d,mobileMemb, fixedDistToCenters=NULL, method="NC", m=2,dnoise=NULL, eta = NULL, alpha=0.001, iter.max=100) {
        METHODS <- c("KM", "FCM", "PCM","NC","HNC" ,"KMdd","NCdd", "HNCdd", "FCMdd", "PCMdd")
        method <- match.arg(method, METHODS)
        if(method=="KM"||method=="KMdd") {
          m=1.0
          dnoise=NULL
          eta=NULL
        }
        else if(method=="FCM"||method=="FCMdd") {
          dnoise=NULL
          eta=NULL
        }
        else if(method=="NC"||method=="NCdd") {
          if(is.null(dnoise)) stop("Must provide a value for dnoise")
          eta = NULL
        }
        else if(method=="HNC"||method=="HNCdd") {
          if(is.null(dnoise)) stop("Must provide a value for dnoise")
          eta = NULL
          m=1.0
        }
        else if(method=="PCM"||method=="PCMdd") {
          if(is.null(eta)) stop("Must provide a vector of values for eta")
          dnoise = NULL
        }
        
        d = as.matrix(d)
        #Number of objects
        n = nrow(d)
        
        #Sets the starting memberships for mobile clusters
        if(is.data.frame(mobileMemb) || is.matrix(mobileMemb)) {
          if(nrow(mobileMemb)!=ncol(d)) {
            stop("The number of rows in mobileMemb must be the same as the number rows and columns of d")
          }		
          u = as.matrix(mobileMemb)
        } else if(is.vector(mobileMemb) && is.numeric(mobileMemb)) {
          u = matrix(0,n,length(mobileMemb))
          for(k in 1:length(mobileMemb)) {
            u[mobileMemb[k],k]=1
          }		
        }
        else {
          stop("Provide a number, a vector of seeds, or membership matrix for mobile clusters")
        }	
        kMov = ncol(u)
        #Sets the fixed cluster memberships
        if(!is.null(fixedDistToCenters)) {
          if(is.data.frame(fixedDistToCenters)) {
            fixedMemb = as.matrix(fixedDistToCenters)
            fixedMemb[] = 0
            kFix = ncol(fixedMemb)
            u = cbind(u,fixedMemb)
          }
          else if(!is.matrix(fixedDistToCenters)) {
            stop("Fixed clusters must be specified as a matrix or a data frame")
          }	
          else {
            fixedMemb = fixedDistToCenters
            fixedMemb[] = 0
            kFix = ncol(fixedMemb)
            u = cbind(u,fixedMemb)
          }
        } else {
          kFix = 0
        }
        #Define vector of medoids
        med = rep(NA,ncol(u))
        
        #Check possibilistic parameters
        if((method=="PCM"||method=="PCMdd") && length(eta)!=(kMov+kFix)) stop("Vector of reference distances (eta) must have a length equal to the number of clusters")
        
        #Add extra (noise) column for NC-related methods
        if(method=="NC"||method=="NCdd"||method=="HNC"||method=="HNCdd") {
          u = cbind(u, vector("numeric",length=n))
        }
        uPrev = matrix(0,nrow=n,ncol=ncol(u))
        
        #Initialize squared distances to fixed centroids
        if(method=="KM"||method=="PCM"||method=="NC"|| method=="HNC"||method=="FCM") {     
          sqdist2cent = matrix(0,nrow=n,ncol=(kMov+kFix))
          if(kFix>0) {
            sqdist2cent[,(kMov+1):(kMov+kFix)] = as.matrix(fixedDistToCenters)^2
          }
        } else { #Initialize distances to fixed medoids
          dist2med = matrix(0,nrow=n,ncol=(kMov+kFix))
          if(kFix>0) {
            dist2med[,(kMov+1):(kMov+kFix)] = as.matrix(fixedDistToCenters)^2
          }
        }
        
        continue = TRUE
        iter = 1
        #iterates until no change in memberships
        while(continue) {
          #1. Update squared distances to mobile centers (centroids for Euclidean-based methods and medoids for the others)
          if(method=="KM"||method=="PCM"||method=="NC"|| method=="HNC"||method=="FCM") {     
            vargeom = vector("numeric", kMov)
            for(k in 1:(kMov)) {
              vargeom[k] = sum((u[,k]^m) %*% (d^2) %*% (u[,k]^m))/(2*sum(u[,k]^m)^2)
              if(is.nan(vargeom[k])) {
                cat(paste0("cluster", k, " vargeom ", vargeom[k],"\n"))
                stop(paste0("NaN vargeom for cluster ", k))
              }
              for(i in 1:n) {
                sqdist2cent[i,k] = (sum((u[,k]^m)*(d[i,]^2))/sum(u[,k]^m))-vargeom[k]
                if(sqdist2cent[i,k]<0) sqdist2cent[i,k]=0
              }
            }
          } else{ 
            for(k in 1:kMov) {
              candidates = 1:n
              excluded = numeric(0)
              # Exclude as candidates objects with membership 1 to other clusters as potential medoids for this cluster
              excluded = which(apply(u[,-k, drop=FALSE], 1 ,max)==1.0)
              if(length(excluded)>0) candidates = candidates[-excluded]
              #Determine medoid
              med[k] = candidates[which.min((u[, k]^m) %*% d[,candidates])]
              dist2med[,k] = d[,med[k]]
            }
          }
          
          #2. compute membership to centroids for mobile and fixed clusters
          if (method=="KM") {
            minC<-apply(sqdist2cent,1,which.min)
            u[,] = 0
            for(k in 1:length(minC)) u[k,minC[k]] = 1.0
          } else if (method=="KMdd") {
            minC<-apply(dist2med,1,which.min)
            u[,] = 0
            for(k in 1:length(minC)) u[k,minC[k]] = 1.0
          } else if(method=="NC") {
            d2cm2<-cbind(sqdist2cent,dnoise^2)
            for(k in 1:ncol(d2cm2)) {
              a<-sweep(d2cm2,1,d2cm2[,k],"/")
              u[,k] = 1/rowSums(a^(-1/(m-1)))
            }
            u[d2cm2==0]=1
          } else if(method=="HNC") {
            d2cm<-cbind(sqdist2cent,dnoise^2)
            u[,] = 0
            minC<-apply(d2cm,1,which.min)
            for(k in 1:length(minC)) {
              u[k,minC[k]] = 1.0
            }
          } else if(method=="HNCdd") {
            d2cm<-cbind(dist2med,dnoise)
            u[,] = 0
            minC<-apply(d2cm,1,which.min)
            for(k in 1:length(minC)) {
              u[k,minC[k]] = 1.0
            }
          } else if(method=="NCdd") {
            d2cm<-cbind(dist2med,dnoise)
            for(k in 1:ncol(d2cm)) {
              a<-sweep(d2cm,1,d2cm[,k],"/")
              u[,k] = 1/rowSums(a^(-1/(m-1)))
            }
            u[d2cm==0]=1
          } else if (method=="FCM") {
            for(k in 1:ncol(sqdist2cent)) {
              a<-sweep(sqdist2cent,1,sqdist2cent[,k],"/")
              u[,k] = 1/rowSums(a^(-1/(m-1)))
            }
            u[sqdist2cent==0]=1
          } else if (method=="FCMdd") {
            d2cm<-dist2med
            for(k in 1:ncol(d2cm)) {
              a<-sweep(d2cm,1,d2cm[,k],"/")
              u[,k] = 1/rowSums(a^(-1/(m-1)))
            }
            u[dist2med ==0]=1
          } else if (method=="PCM") {
            for(k in 1:ncol(sqdist2cent)) u[,k] = 1/(1+((sqdist2cent[,k])/eta[k])^(1/(m-1)))
            u[dist2cent==0]=1
          } else if (method=="PCMdd") {
            for(k in 1:ncol(dist2med)) u[,k] = 1/(1+((dist2med[,k])/eta[k])^(1/(m-1)))
            u[dist2cent==0]=1
          } 	 	
          
          #Check for stopping
          if(iter>2) {
            continue = (max(abs(u-uPrev))>alpha) && (iter<=iter.max) && (max(abs(u-uPrev2))>alpha)
          }   	
          if(continue) {
            uPrev2 = uPrev
            uPrev = u	   	
            iter=iter+1
            if(verbose) cat(".")
          }
        }
        if(method=="FCM" ||  method=="KM") functional = sum(sqdist2cent*(u^m))
        else if(method=="NC"||method=="HNC") functional = sum(sqdist2cent*(u[,-(kMov+kFix+1)]^m))+sum(dnoise^2*u[,kMov+kFix+1]^m)
        else if(method=="FCMdd"||  method=="KMdd") functional = sum(dist2med*(u^m))
        else if(method=="NCdd"||method=="HNCdd") functional = sum(dist2med*(u[,-(kMov+kFix+1)]^m))+sum(dnoise*u[,kMov+kFix+1]^m)
        else if(method=="PCM") {
          functional = 0
          for(k in 1:(kMov+kFix)) functional = functional+sum(sqdist2cent[,k]*(u[,k]^m))+sum(eta[k]*(1-u[,k])^m)
        } else if(method=="PCMdd") {
          functional = 0
          for(k in 1:(kMov+kFix)) functional = functional+sum(dist2med[,k]*(u[,k]^m))+sum(eta[k]*(1-u[,k])^m)
        } 
        if(verbose) cat(paste("\nIterations:", iter,"Functional: ", functional,"\n"))
        
        #Prepare output
        u = as.data.frame(u)   
        if(method=="KM"||method=="PCM"||method=="NC"|| method=="HNC"||method=="FCM") {     
          dist2cent = as.data.frame(sqrt(sqdist2cent))   
        } else {
          dist2cent = as.data.frame(dist2med)   
        }
        for(k in 1:kMov) {
          names(u)[k] = paste("M",k,sep="")
          names(dist2cent)[k] = paste("M",k,sep="")
        }
        if(kFix>1) {
          for(k in (kMov+1):(kMov+kFix)) {
            names(u)[k] = paste("F",k,sep="")
            names(dist2cent)[k] = paste("F",k,sep="")
          }
        }
        if(method=="NC"||method=="NCdd"||method=="HNC"||method=="HNCdd") names(u)[kMov+kFix+1] = "N"
        rownames(u) = rownames(d)
        rownames(dist2cent) = rownames(d)
        size = colSums(u[,1:(kMov+kFix)])
        if(method=="NC"||method=="FCM"||method=="KM"||method=="PCM"||method=="HNC") withinss = colSums((dist2cent^2)*(u[,1:(kMov+kFix)]^m))
        else withinss = colSums((dist2cent)*(u[,1:(kMov+kFix)]^m))
        #Return medoid indices as mobile or fixed centers
        mobileCenters = NULL
        fixedCenters = NULL
        if(method=="KMdd"||method=="FCMdd"||method=="NCdd"||method=="HNCdd"||method=="PCMdd") {
          mobileCenters = med[1:kMov]    
          if(kFix>0) fixedCenters = med[(kMov+1):(kMov+kFix)]
        }
        res = list(mode="dist", method=method, m = m, dnoise = dnoise,
                   eta = eta, memb=u,
                   mobileCenters=mobileCenters, fixedCenters= fixedCenters, 
                   dist2clusters=dist2cent, withinss = withinss, 
                   size=size, functional=functional)
        class(res)<-"vegclust"
        return(res)
      }
    
    x_mat = as.matrix(x)
    n = nrow(x_mat)
    if(is.null(seeds)) seeds = 	1:n
    # Exclude as potential seeds those objects that are at the same position as a previous objects
    toExclude = rep(FALSE, n)
    for(i in 2:n) {
      if(min(x_mat[i,1:(i-1)])==0) toExclude[i] = TRUE
    }
    seeds = seeds[!toExclude]
    #If mobileCenters is a number and nstart>1 perform different random starts
    if(is.vector(mobileMemb) && length(mobileMemb)==1 && is.numeric(mobileMemb)) {
      bestRun = vegclustonedist(x, mobileMemb=sample(seeds, mobileMemb), fixedDistToCenters, method, m,dnoise, eta, alpha, iter.max)
      if(nstart>1) {
        for(i in 2:nstart) {
          run = vegclustonedist(x,mobileMemb=sample(seeds,mobileMemb), fixedDistToCenters, method, m,dnoise,eta, alpha, iter.max)
          if(run$functional<bestRun$functional) {
            bestRun = run
          }
        }
      }		
      return(bestRun)
    } else { #Perform a single run
      return(vegclustonedist(x,mobileMemb, fixedDistToCenters, method, m,dnoise, eta, alpha, iter.max))
    }
  }


