crossmemb<-function(x,y,relativize=TRUE) {
	if(inherits(x, "vegclust") || inherits(x, "vegclass")) x = x$memb
	if(inherits(y, "vegclust") || inherits(y, "vegclass")) y = y$memb
	c=t(x)%*%as.matrix(y)
	if(relativize) c = sweep(c,1,colSums(x),"/")
    return(c)
}