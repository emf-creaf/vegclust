\encoding{UTF-8}
\name{as.memb}
\alias{as.memb}
\title{
Turns into membership matrix
}
\description{
Attempts to turn its cluster vector argument into a membership matrix
}
\usage{
as.memb(cluster)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cluster}{
	A vector indicating the hard membership of each object in \code{x} to a set of groups. Can contain \code{NA} values.
   }
}
\value{
An matrix with as many rows as the length of \code{cluster} and as many columns as different cluster levels. \code{NA} values will have zero membership to all clusters
}
\author{
Miquel De \enc{Cáceres}{Caceres}, CREAF.
}

\seealso{
\code{\link{vegclust}}, \code{\link{vegclass}}
}
\examples{
as.memb(factor(c(1,2,NA)))
}
