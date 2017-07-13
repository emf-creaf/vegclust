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
segmentDistances(d, surveys, distance.type ="directed-segment")
trajectoryDistances(d, surveys, distance.type="DSPD")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
	A site-by-site symmetric matrix or an object of class \code{\link{dist}} containing the distance values between pairs of community states corresponding to the same sites sampled \code{surveys} times. Rows/columns correspond to blocks of surveys: If there are \code{n} sites, then rows/columns are \code{{1 ... n}} for the first survey, \code{{1 ... n}} for the second survey, and so on.
   }
  \item{surveys}{
	Number of sampling surveys.
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
Function \code{trajectoryDistances} returns an object of class \code{\link{dist}} containing the distances between trajectories. Function \code{trajectorySegments} returns a  
}
\author{
Miquel De \enc{Cáceres}{Caceres}, Forest Sciences Center of Catalonia
}
\examples{
##TO BE DONE##
}
