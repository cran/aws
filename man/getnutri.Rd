\name{getnutri}
\alias{getnutri}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Calculate number of grid points in a 3D ellipsoidal neighborhood}
\description{
  This function is for internal use only
}
\usage{
getnutri(radiusq, weights)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{radiusq}{ Square of radius }
  \item{weights}{ axis weights }
}
\details{
}
\value{
  Returns number of grid points (i,j,k) with  $ w[1] i^2 + w[2] j^2 + w[3] k^2 <= radiusq $
}
\references{}
\author{J\"org Polzehl \email{polzehl@wias-berlin.de}}
\note{ ~~further notes~~ }

 ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{awstri}}}

\examples{
}
\keyword{ misc }% at least one, from doc/KEYWORDS
