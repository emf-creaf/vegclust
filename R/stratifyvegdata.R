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
    c1 = cut(treeDataPlot[,size1ColumnId], sizes1)
    colnames(m) = levels(c1)
    c1 = as.numeric(c1)    
    for(i in 1:nrow(treeDataPlot)) {
      if(verbose) cat(paste(i,"\n"))
      isp = which(spcodes==treeDataPlot[i,speciesColumnId])
      if(!cumulative) {
        if(!counts) m[isp,c1[i]] = m[isp,c1[i]]+treeDataPlot[i,abundanceColumnId]
        else m[isp,c1[i]] = m[isp,c1[i]]+1
      } else {
        if(!counts) m[isp,1:c1[i]] = m[isp,1:c1[i]]+treeDataPlot[i,abundanceColumnId]
        else m[isp,1:c1[i]] = m[isp,1:c1[i]]+1  
      }
    }
    return(m)
  }
  doublestratify<-function(treeDataPlot, sizes1, sizes2, spcodes=NULL, speciesColumnId, abundanceColumnId, size1ColumnId, size2ColumnId, cumulative=FALSE, counts=FALSE, verbose=FALSE) {
    if(is.null(spcodes)) spcodes = unique(treeData[,speciesColumnId])
    nsp = length(spcodes)
    nstrata1 = length(sizes1)-1
    nstrata2 = length(sizes2)-1
    m = array(0,dim=c(nsp, nstrata1, nstrata2))
    c1 = cut(treeDataPlot[,size1ColumnId], sizes1)
    c2 = cut(treeDataPlot[,size2ColumnId], sizes2)
    dimnames(m) = list(spcodes, levels(c1), levels(c2))
    c1 = as.numeric(c1)
    c2 = as.numeric(c2)
    for(i in 1:nrow(treeDataPlot)) {
      if(verbose) cat(paste(i,"\n"))
      isp = which(spcodes==treeDataPlot[i,speciesColumnId])
      if(!cumulative){
        if(!counts) m[isp,c1[i],c2[i]] = m[isp,c1[i],c2[i]]+treeDataPlot[i,abundanceColumnId]
        else m[isp,c1[i],c2[i]] = m[isp,c1[i],c2[i]]+1
      } else {
        if(!counts) m[isp,1:c1[i],1:c2[i]] = m[isp,1:c1[i],1:c2[i]]+treeDataPlot[i,abundanceColumnId]
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