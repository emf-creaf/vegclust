\encoding{UTF-8}
\name{trajectories}
\alias{segmentDistances}
\alias{trajectoryDistances}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Distance between directed segments or community trajectories 
}
\description{
Given a distance matrix between community states, functions \code{segmentDistances} and \code{trajectoryDistances} calculate the distance between pairs of directed segments and community trajectories, respectively.
}
\usage{
segmentDistances(d, sites,surveys, distance.type ="directed-segment")
trajectoryDistances(d, sites,surveys, distance.type="DSPD")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
	A site-by-site symmetric matrix or an object of class \code{\link{dist}} containing the distance values between pairs of community states. 
   }
   \item{sites}{
     A vector indicating the site corresponding to each community state.
   }
  \item{surveys}{
     A vector indicating the survey corresponding to each community state (optional).
   }
\item{distance.type}{
The type of distance index to be calculated (De Cáceres et al. in prep). For \code{segmentDistances} the available indices are:
     \itemize{
       \item{\code{Hausdorff}: Hausdorff distance between two segments.}
       \item{\code{directed-segment}: Directed segment distance (default).}
       \item{\code{PPA}: Perpendicular-parallel-angle distance.}
     }
   whereas for \code{trajectoryDistances} the available indices are:
     \itemize{
       \item{\code{Hausdorff}: Hausdorff distance between two trajectories.}
       \item{\code{SPD}: Segment path distance.}
       \item{\code{DSPD}: Directed segment path distance (default).}
     }
  }
}
\value{
Function \code{trajectoryDistances} returns an object of class \code{\link{dist}} containing the distances between trajectories. Function \code{trajectorySegments} returns a list with the following elements:
     \itemize{
       \item{\code{Dseg}: Distance matrix between segments.}
       \item{\code{Dini}: Distance matrix between initial points of segments.}
       \item{\code{Dfin}: Distance matrix between final points of segments.}
       \item{\code{Dinifin}: Distance matrix between initial points of one segment and the final point of the other.}
       \item{\code{Dfinini}: Distance matrix between final points of one segment and the initial point of the other.}
     }

}
\author{
Miquel De \enc{Cáceres}{Caceres}, Forest Sciences Center of Catalonia
}
\examples{
xy<-matrix(0, nrow=6, ncol=2)
xy[2,2]<-1
xy[3,2]<-2
xy[4:6,1] <- 0.5
xy[4:6,2] <- xy[1:3,2]
xy[6,1]<-1
d = dist(xy)
sites = c(1,1,1,2,2,2)
surveys=c(1,2,3,1,2,3)
}
