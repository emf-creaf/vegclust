#' Reshapes community data from individual into stratified form
#' 
#' Function \code{stratifyvegdata} reshapes individual abundance values into species abundance values per size class or combination of size classes. Function \code{as.stratifiedvegdata} checks if the input list has appropriate properties and turns it into an object of class '\code{stratifiedvegdata}'.
#'
#' @param x A data frame containing individual plant data. Individuals are in rows, while measurements are in columns.
#' @param sizes1 A numerical vector containing the breaks for primary size classes in ascending order.
#' @param sizes2 A numerical vector containing the breaks for secondary size classes in ascending order.
#' @param treeSel A logical vector specifying which rows in \code{x} to be used. By default (\code{treeSel = NULL}) all rows are taken.
#' @param spcodes A character vector indicating the codes of species to be used for stratification (species codes beyond those appearing in \code{x} are possible). If \code{spcodes = NULL} then all species in \code{x} are used.
#' @param plotColumn The name of the column in \code{x} that contains plot identifiers.
#' @param speciesColumn The name of the column in \code{x} that contains species names.
#' @param abundanceColumn The name of the column in \code{x} that contains abundance values.
#' @param size1Column The name of the column in \code{x} that contains values for primary size classes.
#' @param size2Column The name of the column in \code{x} that contains values for secondary size classes.
#' @param cumulative A flag to indicate that cumulative abundance profiles or surfaces are desired.
#' @param counts A flag to indicate that the output should be individual counts instead of added abundance values.
#' @param mergeSpecies A flag to indicate that species identity should be ignored. This leads to analyzing the structure of biomass disregarding species identity.
#' @param verbose A logical flag to indicate extra output.
#'
#' @details
#' For each individual (row) in \code{x}, \code{stratifyvegdata} assigns it to the size class (stratum) containing its size. The corresponding abundance value (e.g. crown cover) of the individual is added to the abundance of the corresponding species at the size class (stratum). If \code{sizes2} and \code{size2Column} are supplied, the function assigns each individual (row) in \code{x} to the combination of size classes (e.g. tree height and diameter). 
#' 
#' @returns
#' Both functions return an object of class '\code{stratifiedvegdata}', which is a list of matrices, one for each plot record. Each element (matrix) has as many rows as species and as many columns as size classes (i.e., as many as elements in vector \code{sizes1}). Columns are named starting with 'S' and continuing with the size class (stratum) number. If \code{mergeSpecies=TRUE} then all matrices have a single row (whose name is \code{"all"}). If \code{sizes2} and \code{size2Column} are supplied to \code{stratifyvegdata}, the function returns an object of class '\code{doublestratifiedvegdata}', which is a list of arrays, one for each plot record. Each element (array) has three dimensions corresponding to species, primary sizes (number of elements in in vector \code{sizes1}) and secondary sizes (number of elements in in vector \code{sizes2}). If \code{cumulative=TRUE} then the function returns cumulative abundances (see \code{\link{CAP}} and \code{\link{CAS}}).
#' 
#' @references De \enc{Cáceres}{Caceres}, M., Legendre, P. & He, F. (2013) Dissimilarity measurements and the size structure of ecological communities. Methods in Ecology and Evolution 4: 1167-1177.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF. 
#' 
#' @seealso \code{\link{reshape}}, \code{\link{CAP}}, \code{\link{CAS}}
#' @export
#'
#' @examples
#' ## Load tree data
#' data(treedata)
#' 
#' ## Inspect tree data
#' head(treedata)
#' 
#' ## Define stratum thresholds (4 strata)
#' heights <- seq(0,4, by=0.5)
#' diameters <- seq(0,2, by=0.5)
#' 
#' ## Stratify tree data using heights as structural variable
#' X <- stratifyvegdata(treedata, sizes1=heights, plotColumn="plotID",
#'                      speciesColumn="species", size1Column="height", counts=TRUE)
#' 
#' ## Inspect the second plot record
#' X[[2]]
#' 
#' ## Stratify tree data using heights as structural variable and cover as abundance
#' Y <- stratifyvegdata(treedata, sizes1=heights, plotColumn="plotID",
#'                      speciesColumn="species", size1Column="height", 
#'                      abundanceColumn="cover")
#' Y[[2]]
#' 
#' ## Stratify tree data using heights and diameters as structural variables
#' Z <- stratifyvegdata(treedata, sizes1=heights, sizes2=diameters, plotColumn="plotID",
#'                      speciesColumn="species", size1Column="height", size2Column="diam",
#'                      counts=TRUE)
#' Z[[2]]
#' 
stratifyvegdata<-function(x,sizes1, sizes2 = NULL, treeSel=NULL, spcodes=NULL, plotColumn="plot", 
                          speciesColumn = "species", abundanceColumn="abundance", size1Column = "size", 
                          size2Column = NULL, cumulative = FALSE, counts=FALSE, mergeSpecies=FALSE, verbose=FALSE) {
  treeData = as.data.frame(x)
  plotColumnId = which(names(treeData)==plotColumn)
  abundanceColumnId = which(names(treeData)==abundanceColumn)
  size1ColumnId = which(names(treeData)==size1Column)
  doublestratify = FALSE
  if(!is.null(size2Column) || (!is.null(sizes2))) {
    doublestratify = TRUE
    if(is.null(sizes2)) stop("sizes2 must be specified for double stratification")
    if(is.null(size2Column)) stop("size2Column must be specified for double stratification")
    size2ColumnId = which(names(treeData)==size2Column)
  }
  speciesColumnId = which(names(treeData)==speciesColumn)
  if(mergeSpecies) treeData[,speciesColumnId] = "allspecies"
  if(is.null(spcodes)) spcodes =sort(unique(treeData[,speciesColumnId]))
  if(!is.null(treeSel)) {
    treeData = treeData[treeSel,]
  }
  stratify<-function(treeDataPlot, sizes, spcodes=NULL, speciesColumnId, abundanceColumnId, sizeColumnId, cumulative=FALSE, counts=FALSE, verbose=FALSE) {
    if(is.null(spcodes)) spcodes = unique(treeData[,speciesColumnId])
    nsp = length(spcodes)
    nstrata = length(sizes)-1
    m = matrix(0,nrow=nsp, ncol=nstrata)
    rownames(m) = spcodes
    c1 = cut(treeDataPlot[,sizeColumnId], sizes, include.lowest = TRUE)
    if(sum(is.na(c1))>0) stop("Some values are not included within size classes. Revise size class definition")
    colnames(m) = levels(c1)
    c1 = as.numeric(c1)    
    for(i in 1:nrow(treeDataPlot)) {
      isp = which(spcodes==treeDataPlot[i,speciesColumnId])
      if(verbose) cat(paste(i,"_",isp,"_",c1[i],"_",m[isp, c1[i]], ": ", treeDataPlot[i,abundanceColumnId],"\n"))
      if(!cumulative) {
        if(!counts) m[isp, c1[i]] = m[isp,c1[i]]+ as.numeric(treeDataPlot[i,abundanceColumnId])
        else m[isp,c1[i]] = m[isp,c1[i]]+1
      } else {
        if(!counts) m[isp,1:c1[i]] = m[isp,1:c1[i]]+as.numeric(treeDataPlot[i,abundanceColumnId])
        else m[isp,1:c1[i]] = m[isp,1:c1[i]]+1  
      }
    }
    return(m)
  }
  doublestratify<-function(treeDataPlot, sizes1, sizes2, spcodes=NULL, speciesColumnId, abundanceColumnId, 
                           size1ColumnId, size2ColumnId, cumulative=FALSE, counts=FALSE, verbose=FALSE) {
    if(is.null(spcodes)) spcodes = unique(treeData[,speciesColumnId])
    nsp = length(spcodes)
    nstrata1 = length(sizes1)-1
    nstrata2 = length(sizes2)-1
    m = array(0,dim=c(nsp, nstrata1, nstrata2))
    c1 = cut(treeDataPlot[,size1ColumnId], sizes1, include.lowest = TRUE)
    if(sum(is.na(c1))>0) stop("Some values are not included within size1 classes. Revise size1 class definition")
    c2 = cut(treeDataPlot[,size2ColumnId], sizes2, include.lowest = TRUE)
    if(sum(is.na(c2))>0) stop("Some values are not included within size2 classes. Revise size2 class definition")
    dimnames(m) = list(spcodes, levels(c1), levels(c2))
    c1 = as.numeric(c1)
    c2 = as.numeric(c2)
    for(i in 1:nrow(treeDataPlot)) {
      if(verbose) cat(paste(i,"\n"))
      isp = which(spcodes==treeDataPlot[i,speciesColumnId])
      if(!cumulative){
        if(!counts) m[isp,c1[i],c2[i]] = m[isp,c1[i],c2[i]]+as.numeric(treeDataPlot[i,abundanceColumnId])
        else m[isp,c1[i],c2[i]] = m[isp,c1[i],c2[i]]+1
      } else {
        if(!counts) m[isp,1:c1[i],1:c2[i]] = m[isp,1:c1[i],1:c2[i]]+as.numeric(treeDataPlot[i,abundanceColumnId])
        else m[isp,1:c1[i],1:c2[i]] = m[isp,1:c1[i],1:c2[i]]+1
      }
    }
    return(m)
  }
  
  if(!is.null(size2Column)){ #double stratification
    X = lapply(split(treeData,treeData[,plotColumnId]),
               FUN=doublestratify,sizes1=sizes1, sizes2 =sizes2,
               spcodes=spcodes, speciesColumnId = speciesColumnId, abundanceColumnId =abundanceColumnId, size1ColumnId=size1ColumnId,size2ColumnId=size2ColumnId, 
               cumulative=cumulative, counts=counts, verbose=verbose)
    if(!cumulative) class(X)<-c("doublestratifiedvegdata","list")
    else class(X)<-c("CAS","list")  
  } else {
    X = lapply(split(treeData,treeData[,plotColumnId]),
               FUN=stratify,sizes=sizes1, 
               spcodes=spcodes, speciesColumnId = speciesColumnId, abundanceColumnId =abundanceColumnId, sizeColumnId=size1ColumnId, 
               cumulative = cumulative, counts=counts, verbose=verbose)
    if(!cumulative) class(X)<-c("stratifiedvegdata","list")
    else class(X)<-c("CAP","list")
  }
  return(X)
}



#' @rdname stratifyvegdata
#' @export
#' @param X A list with as many elements as plot records. Each element should be of class 'matrix' or 'data.frame' with species in rows and strata in columns. Furthermore, the number of rows (species) and columns (strata) should be the same for all elements.
as.stratifiedvegdata<-function(X) {
  if(!inherits(X,"list")) stop("Input should be of class 'list'")
  n = length(X)
  if(n==0) stop("The list should be of positive length")
  for(i in 1:n) if(!inherits(X[[i]],"matrix") && !inherits(X[[i]],"data.frame")) stop("The list should contain 'data.frame' or 'matrix' objects only")
  p = nrow(X[[1]])
  s = ncol(X[[1]])
  for(i in 2:n) {
    if(nrow(X[[i]])!=p) stop("All elements should have the same number of rows (species)")
    if(ncol(X[[i]])!=s) stop("All elements should have the same number of columns (strata)")
  }  
  class(X)<-c("stratifiedvegdata","list")
  if(is.null(names(X))) names(X)<-1:length(X)
  return(X)
}