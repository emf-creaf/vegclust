#' Constancy table of a classification
#' 
#' Allows studying the constancy table (i.e. the frequency of species in each class) of a classification represented in the form of a membership data matrix.
#'
#' @param x Community data, a site by species data frame.
#' @param memb An site-by-group matrix indicating the (hard or fuzzy) membership of each object in \code{x} to a set of groups.
#'
#' @returns
#' Function \code{clustconst} returns an object of type 'clustconst', in fact a data frame with the constancy value of each species (rows) on each cluster (column).
#' 
#' @details
#' The constancy value of a species in a vegetation unit is the relative frequency of occurrence of the species in plot records that belong to the unit. In case of a fuzzy vegetation unit the constancy value is the sum of memberships of sites that contain the species divided by the sum of memberships of all sites. Use the 'summary' function to obtain information about: (1) which species are more frequent on a given vegetation unit; (2) which vegetation units have higher frequencies of a given target species. Additionally, the 'summary' function can sort a constancy table if \code{mode="all"} and \code{sort=TRUE} are indicated.
#' 
#' @author 
#' Miquel De \enc{Cáceres}{Caceres}, CREAF
#' 
#' @seealso 
#' \code{\link{vegclust}}, \code{\link{kmeans}}
#' 
#' @export
#' 
#' @name clustconst
#'
#' @examples
#' ## Loads stats
#' library(stats)  
#' 
#' ## Loads data  
#' data(wetland)
#' 
#' ## This equals the chord transformation 
#' ## (see also \code{\link{decostand}} in package 'vegan')
#' wetland.chord <- as.data.frame(sweep(as.matrix(wetland), 1, 
#'                                      sqrt(rowSums(as.matrix(wetland)^2)), "/"))
#' 
#' ## Performs a K-means clustering
#' wetland.km <- kmeans(wetland.chord, centers=3, nstart=10)
#' 
#' ## Gets constancy table of KM (i.e. hard) clusters
#' c <- clustconst(wetland.chord, memb=as.memb(wetland.km$cluster))
#' 
#' ## Prints constancy values ordered and store the result in d
#' d <- summary(c, mode="all")
#' 
#' ## Prints the most frequent species in the first cluster
#' summary(c, mode="cluster", name=names(c)[1])
#' 
clustconst <- function(x,memb) {
	pa<-ifelse(x>0,1,0)
	v<-sweep(t(pa)%*%as.matrix(memb),2,colSums(memb),"/")
	v = as.data.frame(v)
	if(is.data.frame(memb)) names(v) = names(memb)
	if(is.data.frame(x)) row.names(v)=names(x)
	class(v)<-c("clustconst","data.frame")
	return(v)		
}

#' @rdname clustconst
#' @param object An object of class 'clustconst'.
#' @param mode Use \code{mode="all"} to print the constancy table, \code{mode="cluster"} to print constancy values for one cluster, and \code{mode="species"}, to print constancy values for one species.
#' @param name A string with the name of a cluster (in \code{mode="cluster"}), or the name of a species (in \code{mode="species"}).
#' @param sort A flag to indicate whether constancy table should be sorted in descending order.
#' @param minconst A threshold used to limit the values shown.
#' @param digits The number of digits for rounding.
#' @param ... Additional parameters for summary (actually not used).
#' 
#' @export
summary.clustconst<-function(object, mode="all", name=NULL, sort=TRUE, minconst=0.5, digits=3,...) {
  x <- as.data.frame(object)
  if(mode=="cluster") {
    clustIndex = which(names(x)==name)
    spnames = row.names(x)
    values = x[,clustIndex]
    if(sort) {
      o = order(x[, clustIndex],decreasing=TRUE)
      spnames = spnames[o]
      values = values[o]
    }
    spIndices = which(values>= minconst)
    for(i in 1:length(spIndices)) {
      cat(paste(spnames[spIndices[i]],format(round(values[spIndices[i]],digits=digits),nsmall=digits),"\n"))
    }
  }
  else if(mode=="species") {
    spIndex = which(row.names(x)==name)
    clnames = names(x)
    values = x[spIndex,]
    if(sort) {
      o = order(x[spIndex,],decreasing=TRUE)
      clnames = clnames[o]
      values = values[o]
    }	
    clIndices = which(values>=minconst)
    for(i in 1:length(clIndices)) {
      cat(paste(clnames[clIndices[i]],format(round(values[clIndices[i]],digits=digits),nsmall=digits),"\n"))
    }
  }
  else if(mode=="all") {
    x =x[apply(x,1,max)>minconst,]
    if(sort) {
      y = sweep(x,1,rowSums(x),"/")
      oc = order(colSums(y), decreasing=TRUE)
      y = y[,oc]
      x = x[,oc]
      o = integer(nrow(x))
      t = 0
      used=logical(nrow(x))
      for(k in 1:ncol(x)) {
        indices = which(apply(y,1,which.max)==k & !used)
        if(length(indices)>0) {
          o[(t+1):(t+length(indices))] = indices[order(x[indices,k], decreasing=TRUE)]
          t = t+length(indices)
          
          cat(paste("------------",names(x)[k],"-------------\n"))
          print(format(round(x[indices[order(x[indices,k], decreasing=TRUE)],],digits=digits), nsmall=digits))
          used[indices] = TRUE
        }
      }
      if(sum(!used)>0) {
        o[(t+1):(t+sum(!used))] = which(!used)
        cat(paste("------------ REMAINING -------------\n"))
        print(format(round(x[which(!used),],digits=digits), nsmall=digits))
      }
      x = x[o,]
    }
    invisible(x)
  }
}
