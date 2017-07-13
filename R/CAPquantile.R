CAPquantile<-function(CAP, q = 0.5, y = NULL) {
  quantileCAP<-function(x, q = 0.5) {
    res = x[[1]]
    v = numeric(length(x))
    for(r in 1:nrow(res)) {
      for(c in 1:ncol(res)) {
        for(i in 1:length(x)) v[i] = x[[i]][r,c]
        res[r,c] = quantile(v, probs = q)
      }
    }
    return(res)
  }
  if(!is.null(y)) {
    if(is.vector(y) || is.factor(y)) {
      CC = lapply(split(CAP, as.factor(y)), FUN = quantileCAP, q)
    }    
  } else {
    CC = list(qCAP=quantileCAP(CAP, q))
  }
  if(inherits(CAP, "CAP")) class(CC)<-c("CAP","list")
  else if(inherits(CAP, "stratifiedvegdata")) class(CC)<-c("stratifiedvegdata","list")
  return(CC)  
}