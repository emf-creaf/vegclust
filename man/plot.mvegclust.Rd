\encoding{UTF-8}
\name{plot.mvegclust}
\alias{plot.mvegclust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots clustering results 
}
\description{
Create plots used to study vegclust clustering results for an increasing number of clusters
}
\usage{
\method{plot}{mvegclust}(x, type="hnc", excludeFixed=TRUE, verbose=FALSE, ylim=NULL, 
                         xlab=NULL, ylab=NULL, maxvar=0.6, minsize=20,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object returned from functions \code{\link{hier.vegclust}} or \code{\link{random.vegclust}}.}
  \item{type}{A string indicating the type of plot desired. Current accepted values are "hnc","hmemb","var","hcs" and "valid".}
  \item{excludeFixed}{A flag to indicate whether clusters with fixed centroids should be excluded from plots.}
  \item{verbose}{A flag to print extra information.}
  \item{ylim}{A vector with the limits for the y axis.}
  \item{xlab}{String label for the x axis.}
  \item{ylab}{String label for the y axis.}
  \item{maxvar}{Maximum cluster variance allowed for the \code{type="valid"} plot.}
  \item{minsize}{Minimum cluster size allowed for the \code{type="valid"} plot.}
  \item{...}{Additional plotting parameters.}
}
\value{
Different information is returned depending on the type of plot chosen.
}
\author{
Miquel De \enc{Cáceres}{Caceres}, Forest Science Center of Catalonia
}
\examples{
## Loads data  
data(wetland)
  
## This equals the chord transformation 
## (see also \code{\link{decostand}} in package 'vegan')
wetland.chord = as.data.frame(sweep(as.matrix(wetland), 1, 
                              sqrt(rowSums(as.matrix(wetland)^2)), "/"))

## Create noise clustering from hierarchical clustering at different number of clusters
wetland.hc = hclust(dist(wetland.chord),method="ward") 
wetland.nc = hier.vegclust(wetland.chord, wetland.hc, cmin=2, cmax=5, m = 1.2, 
                           dnoise=0.75, method="NC")

## Plot changes in the number of objects falling into the noise cluster
plot(wetland.nc, type="hnc")

## Plots the number of objects falling into "true" clusters, 
## the number of objects considered intermediate, 
## and the number of objects falling into the noise
plot(wetland.nc, type="hmemb")

## Plot minimum, maximum and average cluster size
plot(wetland.nc, type="hcs")

## Plot minimum, maximum and average cluster variance
plot(wetland.nc, type="var")

## Plot number of groups with high variance, low membership or both
plot(wetland.nc, type="valid")
}
