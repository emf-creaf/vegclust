\encoding{UTF-8}
\name{concordance}
\alias{concordance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Concordance between two classifications
}
\description{
Computes an index to compare two classifications.
}
\usage{
concordance(x, y, method="adjustedRand",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,y}{ Classification vector or membership matrix. Alternatively, objects of type \code{\link{vegclust}} or \code{\link{vegclass}}.}
  \item{method}{A string vector to indicate the desired indices (see details).}
  \item{...}{
    Additional parameters for function \code{\link{defuzzify}}, which will be called if \code{x} or \code{y} are of type \code{matrix}, \code{\link{vegclust}} or \code{\link{vegclass}}.
  }
}
\details{
Several indices for comparison of partitions are available:
\itemize{
  \item{\code{method="Rand"}: Rand (1971) index.}
  \item{\code{method="adjustedRand"}: Rand index adjusted for random effects (Hubert & Arabie 1985).}
  \item{\code{method="Wallace"}: Wallace (1983) index (for asymmetrical comparisons). This index (and its adjusted version) is useful to quantify how much \code{x} is nested into \code{y}.}
  \item{\code{method="adjustedWallace"}: Wallace index adjusted for random effects (Pinto et al. 2008).}
}
}
\references{
Hubert, L. & Arabie, P. (1985). Comparing partitions. Journal of Classification, 2, 193–218.

Pinto, F.R., Melo-Cristino, J. & Ramirez, M. (2008). A confidence interval for the wallace coefficient of concordance and its application to microbial typing methods. PLoS ONE, 3.

Rand, W.M. (1971). Objective Criteria for the Evaluation of Clustering Methods. Journal of the American Statistical Association, 66, 846–850.

Wallace, D.L. (1983). A method for comparing two hierarchical clusterings: Comment. Journal of the American Statistical Association, 78, 569–576.
}
\value{
  A numeric vector with the desired index values.
}
\author{
Miquel De \enc{Cáceres}{Caceres}, Forest Science Center of Catalonia
}

\seealso{
\code{\link{vegclust}}, \code{\link{vegclass}}, \code{\link{defuzzify}}, \code{\link[mclust]{adjustedRandIndex}}
}
%\examples{
%}
