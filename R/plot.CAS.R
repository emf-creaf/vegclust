plot.CAS<-function(x, plot=NULL, species=NULL, sizes1=NULL, sizes2 = NULL, palette = colorRampPalette( c( "light blue","light green","white", "yellow","orange","red") ), zlim=NULL, ...){
  if(!inherits(x,"CAS")) stop("Input should be of class 'CAS'")
  if(is.null(plot)) plot=names(x)[1]
  if(is.null(species)) species = dimnames(x[[1]])[[1]][1]
  a<-x[[plot]] #Select plot to draw
  z <-a[species,,]
  if(is.null(sizes1)) sizes1 = 1:nrow(z)
  if(is.null(sizes2)) sizes2 = 1:ncol(z)  
  # Generate the desired number of colors from this palette
  nbcol <- 100
  color <- palette(nbcol)
  # Compute the z-value at the facet centres
  nrz <- nrow(z)
  ncz <- ncol(z)
  zfacet <- (z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz])/4
  # Recode facet z-values into color indices
  if(!is.null(zlim)) facetcol <- cut(zfacet, 
                                     breaks=seq(zlim[1], zlim[2], length.out=(nbcol+1)),
                                     include.lowest=TRUE)
  else facetcol <- cut(zfacet, nbcol)  
  if(!is.null(zlim)) persp(x=sizes1, y=sizes2,z, col = color[facetcol], ticktype="detailed", 
        theta=130, zlim=zlim, ...)
  else persp(x=sizes1, y=sizes2,z, col = color[facetcol], ticktype="detailed", 
             theta=130, ...)
}