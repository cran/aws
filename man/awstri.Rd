\name{awstri}
\alias{awstri}
\title{Tri-dimensional Adaptive Weights Smoothing}
\description{Performes tri dimensional Adaptive Weigths Smoothing}
\usage{awsbi(y, lambda=3, gamma=1.3, eta =4, s2hat = NULL, kstar = 
length(radii),
       rmax=max(radii), weight = c(1,1,1), radii = 
c((1:4)/2,2.3,(5:12)/2,7:9,10.5,12,13.5),
       control="dyadic"}
\arguments{
\item{y}{array of observed values}
\item{lambda}{main smoothing parameter (should be approximately 3)}
\item{gamma}{allow for increase of variances during iteration by factor gamma 
(!! gamma >=1) }
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
\references{Polzehl, J. and Spokoiny, S. (2000). Adaptive Weights Smoothing 
with applications
to image restoration. J.R.Statist.Soc. B, 62, Part 2, pp.335-354}
\note{The function assumes that the data are given on a 3D-grid 
corresponding to the dimensionality of y}
\author{J\"org Polzehl \email{polzehl@wias-berlin.de}}
\seealso{\code{\link{awsbi}},\code{\link{awsuni}}}

\examples{
xy_rbind(rep(0:30,31),rep(0:30,rep(31,31)))
w3_array(0,c(31,31,31))
w3[4:28,4:28,4:28]_1
dim(w3)_c(961,31)
w3[((xy[1,]-15)^2+(xy[2,]-15)^2)<=144,16]_0
for(i in 1:12) {
   r2_144-i*i
   w3[((xy[1,]-15)^2+(xy[2,]-15)^2)<=r2,16+c(-i,i)]_0
}
dim(w3)_c(31,31,31)
w3[10:22,10:22,10:22]_1
dim(w3)_c(961,31)
w3[((xy[1,]-15)^2+(xy[2,]-15)^2)<=36,16]_0
for(i in 1:6) {
   r2_36-i*i
   w3[((xy[1,]-15)^2+(xy[2,]-15)^2)<=r2,16+c(-i,i)]_0
}
dim(w3)_c(31,31,31)
sigma_.4
y_w3+rnorm(w3,0,sigma)
yhat_awstri(y,rmax=10)
par(mfrow=c(1,3))
for(i in 1:31){
image(y[,,i],col=gray((0:255)/255))
title(paste("Layer",i,"  Noisy cube"))
image(yhat$yhat[,,i],zlim=range(y),col=gray((0:255)/255))
title(paste("Layer",i,"  AWS reconstruction"))
image(w3[,,i],zlim=range(y),col=gray((0:255)/255))
title(paste("Layer",i,"  Original Cube"))
#readline()
}
rm(y,yhat,w3,xy)
}
\keyword{regression}
\keyword{nonparametric}
\keyword{smooth}

