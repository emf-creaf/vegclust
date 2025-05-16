#' Turns into membership matrix
#'
#' Attempts to turn its cluster vector argument into a membership matrix
#' @param cluster A vector indicating the hard membership of each object in \code{x} to a set of groups. Can contain \code{NA} values.
#'
#' @returns
#' An matrix with as many rows as the length of \code{cluster} and as many columns as different cluster levels. \code{NA} values will have zero membership to all clusters
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF.
#' @seealso \code{\link{vegclust}}, \code{\link{vegclass}}
#' 
#' @export
#'
#' @examples
#' as.memb(factor(c(1,2,NA)))
as.memb<-function(cluster){
   cln =levels(as.factor(cluster))
   k = length(cln)
   u = matrix(0,nrow=length(cluster),ncol=k)    
   for(i in 1:k) {
   	u[,i]=ifelse(cluster==cln[i],1,0)
   } 
   u[is.na(u)]<-0
   if(!is.null(names(cluster))) rownames(u) = names(cluster)
   colnames(u) = cln
   return(u)
}