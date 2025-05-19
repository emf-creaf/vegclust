#' Draws a cummulative abundance surface
#' 
#' Create plots used to inspect one or more cumulative abundance surfaces.
#'
#' @param x An object of class \code{\link{CAS}}.
#' @param plot A string indicating the plot record whose surface is to be drawn.
#' @param species A string indicating the species whose profile is to be drawn.
#' @param sizes1 A vector containing the size values associated to each primary size class. If \code{NULL} the x-axis will be defined using the primary size class order in \code{x}.
#' @param sizes2 A vector containing the size values associated to each secondary size class. If \code{NULL} the y-axis will be defined using the secondary size class order in \code{x}.
#' @param palette Color palette for z values.
#' @param zlim The limits for the z-axis.
#' @param ... Additional plotting parameters for function \code{\link{persp}}.
#'
#' @references
#' De \enc{Cáceres}{Caceres}, M., Legendre, P. & He, F. (2013) Dissimilarity measurements and the size structure of ecological communities. Methods in Ecology and Evolution 4: 1167-1177.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' 
#' @seealso \code{\link{CAS}}, \code{\link{persp}}
#' 
#' @export
#'
#' @examples
#' ## Create synthetic tree data
#' pl <- rep(1,100) # All trees in the same plot
#' sp <- ifelse(runif(100)>0.5,1,2) # Random species identity (species 1 or 2)
#' h <- rgamma(100,10,2) # Heights (m)
#' d <- rpois(100, lambda=h^2) # Diameters (cm)
#' m <- data.frame(plot=pl,species=sp, height=h,diameter=d) 
#' m$ba <- pi*(m$diameter/200)^2
#' print(head(m))
#' 
#' ## Size classes
#' heights <- seq(0,4, by=.25)^2 # Quadratic classes
#' diams <- seq(0,130, by=5) # Linear classes
#' 
#' ## Stratify tree data
#' X <- stratifyvegdata(m, sizes1=heights, sizes2=diams, 
#'                      plotColumn = "plot", speciesColumn = "species", 
#'                      size1Column = "height", size2Column = "diameter", 
#'                      abundanceColumn = "ba")
#' 
#' ## Build cummulative abundance surface
#' Y <- CAS(X)
#' 
#' ## Plot the surface of species '1' in plot '1' using heights and diameters
#' plot(Y, species=1, sizes1=heights[-1], xlab="height (m)", 
#'      ylab="diameter (cm)", sizes2=diams[-1], zlab="Basal area (m2)",
#'      zlim = c(0,6), main="Species 1")
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