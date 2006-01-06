\name{awstri}
\alias{awstri}
\title{Three-dimensional Adaptive Weights Smoothing}
\description{Performes three-dimensional Adaptive Weigths Smoothing (depreciated version, use aws instead)}
\usage{awstri(y, lambda=3, gamma=1.3, eta =4, s2hat = NULL, kstar = length(radii),
       rmax=max(radii), weight = c(1,1,1), radii = c((1:4)/2,2.3,(5:12)/2,7:9,10.5,12,13.5), control="dyadic")}
\arguments{
\item{y}{array of observed values}
\item{lambda}{main smoothing parameter (should be approximately 3)}
\item{gamma}{allow for increase of variances during iteration by factor gamma 
 (!!\code{gamma >=1}) }
\item{eta}{main control parameter (should be approximately 4)   }
\item{s2hat}{initial variance estimate (if available,
           can be either a number (homogeneous case), a matrix of same 
dimension 
           as y (inhomogeneous variance) or NULL (a homogeneous variance 
estimate
           will be generated in this case)}
\item{kstar}{maximal number of iterations to perform, actual number may be 
smaller depending
             on parameters radii and rmax}
\item{weight}{weights used for distances, determining elliptical neighborhoods}
\item{radii}{radii of circular neighbourhoods used}
\item{rmax}{maximal radius of neighborhood to be used, may change kstar}
\item{control}{the control step is performed in either a dyadic sceme
           ("dyadic") or using all previous estimates (otherwise)}
}

\value{
A list with components
\item{yhat}{estimates of the regression function (matrix corresponding to the 
y's)}
\item{shat}{estimated standard deviations of yhat (conditional on the chosen 
weights)}
\item{args}{Main arguments supplied to awstri}
}
\references{Polzehl, J. and Spokoiny, V. (2000). \emph{Adaptive Weights Smoothing 
with applications
to image restoration}, J.R.Statist.Soc. B, 62, Part 2, pp.335-354}
\note{The function assumes that the data are given on a 3D-grid 
corresponding to the dimensionality of y. 
This function is superseded by function aws and will be 
removed in the next mayor version of the package.}
\author{Joerg Polzehl \email{polzehl@wias-berlin.de}}
\seealso{\code{\link{aws}}, \code{\link{awsbi}},\code{\link{awsuni}}}

\examples{
xy <- rbind(rep(0:30,31),rep(0:30,rep(31,31)))
w3 <- array(0,c(31,31,31))
w3[4:28,4:28,4:28] <- 1
dim(w3) <- c(961,31)
w3[((xy[1,]-15)^2+(xy[2,]-15)^2)<=144,16] <- 0
for(i in 1:12) {
   r2 <- 144-i*i
   w3[((xy[1,]-15)^2+(xy[2,]-15)^2)<=r2,16+c(-i,i)] <- 0
}
dim(w3) <- c(31,31,31)
w3[10:22,10:22,10:22] <- 1
dim(w3) <- c(961,31)
w3[((xy[1,]-15)^2+(xy[2,]-15)^2)<=36,16] <- 0
for(i in 1:6) {
   r2 <- 36-i*i
   w3[((xy[1,]-15)^2+(xy[2,]-15)^2)<=r2,16+c(-i,i)] <- 0
}
dim(w3) <- c(31,31,31)
sigma <- .4
y <- w3+rnorm(w3,0,sigma)
#  increase rmax for better results
yhat <- awstri(y,rmax=2)
rm(y,yhat,w3,xy)
}
\keyword{regression}
\keyword{nonparametric}
\keyword{smooth}
