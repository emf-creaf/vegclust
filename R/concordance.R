concordance<-function(x,y, method="adjustedRand",...){
  method = match.arg(method,choices=c("Rand","adjustedRand","Wallace","adjustedWallace"), 
                      several.ok=TRUE)
  if(inherits(x, "vegclust") || inherits(x, "vegclass") || inherits(x, "matrix")) x<-defuzzify(x,...)$cluster
  if(inherits(y, "vegclust") || inherits(y, "vegclass") || inherits(y, "matrix")) y<-defuzzify(y,...)$cluster
  if (length(x) != length(y)) 
    stop("arguments must be classifications of the same number of objects")
  tab <- table(x, y)
  N<-length(x)
  if (all(dim(tab) == c(1, 1))) return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  I<-numeric(length(method))
  for(i in 1:length(method)) {
    if(method[i]=="Rand") {
      I[i]<-(a+d)/(a+b+c+d)
    } else if(method[i]=="adjustedRand"){
      I[i] <- (a - (a + b) * (a + c)/(a + b + c + d))/
        ((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
    } else if(method[i]=="Wallace") {
      I[i] <- a/(a+b)
    } else if(method[i]=="adjustedWallace"){
      W <-a/(a+b)
      nj = colSums(tab)
      SID = (sum(nj*(nj-1))/(N*(N-1)))
      I[i] = (W-SID)/(1-SID)
    }
  }
  names(I)<-method
  return(I)
}