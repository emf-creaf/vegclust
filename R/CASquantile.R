CASquantile<-function(CAS, q = 0.5, y = NULL) {
  quantileCAS<-function(x, q = 0.5) {
    res = x[[1]]
    v = numeric(length(x))
    d = dim(res)
    for(r in 1:d[1]) {
      for(c in 1:d[2]) {
        for(s in 1:d[3]) {
          for(i in 1:length(x)) v[i] = x[[i]][r,c,s]
          res[r,c,s] = quantile(v, probs = q)
        }
      }
    }
    return(res)
  }
  if(!is.null(y)) {
    if(is.vector(y) || is.factor(y)) {
      CC = lapply(split(CAS, as.factor(y)), FUN = quantileCAS, q)
    }    
  } else {
    CC = list(qCAS=quantileCAS(CAS, q))
  }
  if(inherits(CAS, "CAS")) class(CC)<-c("CAS","list")
  return(CC)  
}