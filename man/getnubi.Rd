\name{getnubi}
\alias{getnubi}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Calculate number of grid points in a bivariate elliptical neighborhood}
\description{
  This function is for internal use only
}
\usage{
getnubi(radiusq, weights)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{radiusq}{ Square of radius }
  \item{weights}{ axis weights }
}
\details{
}
\value{
  Returns number of grid points (i,j) with  $ w[1] i^2 + w[2] j^2 <= radiusq $
}
\author{J\"org Polzehl \email{polzehl@wias-berlin.de}}
\seealso{\code{\link{awstri}}}
\examples{
}
\keyword{misc}% at least one, from doc/KEYWORDS
