\name{awsbi}
\alias{awsbi}
\title{Two-dimensional Adaptive Weights Smoothing}
\description{Performes two dimensional Adaptive Weigths Smoothing}
\usage{awsbi(y, lambda=3, gamma=1.3, eta =4, s2hat = NULL, kstar = 
length(radii),
              rmax=max(radii), radii = c((1:8)/2,4.4,5.,(6:10),(6:10)*2),
              graph = F,u0 = NULL, control="dyadic", demomode=F, 
colors=gray((0:255)/255))}
\arguments{
\item{y}{matrix of observed values}
\item{lambda}{main smoothing parameter (should be approximately 3)}
\item{gamma}{allow for increase of variances during iteration by factor gamma 
(!! gamma >=1) }
\item{eta}{main control parameter (should be approximately 4)   }
\item{s2hat}{initial variance estimate (if available,
           can be either a number (homogeneous case), a matrix of same
           dimension as y (inhomogeneous variance) or NULL (a homogeneous 
           variance estimate will be generated in this case)}
\item{kstar}{maximal number of iterations to perform, actual number may be 
           smaller depending on parameters radii, rmax and eps}
\item{radii}{radii of circular neighbourhoods used}
\item{rmax}{maximal radius of neighborhood to be used, may change kstar}
\item{graph}{logical, if T progress (for each iteration) is illustrated 
           grahically, if F the program runs until the final estimate is 
           obtained (much faster !!!)}
\item{u0}{allows for submission of "true" values for illustration and test 
           purposes; only if graph=T,  MSE and MAE are reported for each 
           iteration step}
\item{control}{the control step is performed in either a dyadic sceme
           ("dyadic") or using all previous estimates (otherwise)}
\item{demomode}{if T the function will wait for user input after each 
           iteration; only if graph=T }
\item{colors}{color sceme to be used for images}
}
\value{A list with components
\item{yhat}{estimates of the regression function (matrix corresponding to the 
y's)}
\item{shat}{estimated standard deviations of yhat (conditional on the chosen 
weights)}
\item{nu}{maximal number of design points in neighborhood used}
\item{args}{Main arguments supplied to awsbi}
}
\references{Polzehl, J. and Spokoiny, S. (2000). Adaptive Weights Smoothing 
with applications
to image restoration. J.R.Statist.Soc. B, 62, Part 2, pp.335-354}
\note{The function assumes that the data are given on a 2D-grid 
corresponding to the dimensionality of y}
\author{J\"org Polzehl \email{polzehl@wias-berlin.de}}
\seealso{\code{\link{awsuni}},\code{\link{awstri}}}

\examples{
xy_rbind(rep(0:255,256),rep(0:255,rep(256,256)))
indw_c(1:12,29:48,73:100,133:168,209:256)
w0_matrix(rep(F,256*256),ncol=256)
w0[indw,]_T
w0[,indw]_!w0[,indw]
w0_w0-.5
w0[((xy[1,]-129)^2+(xy[2,]-129)^2)<=10000&((xy[1,]-129)^2+(xy[2,]-129)^2)>=4900]_0
w0[abs(xy[1,]-xy[2,])<=20&((xy[1,]-129)^2+(xy[2,]-129)^2)<4900]_0
w0[((xy[1,]-225)^2+2*(xy[2,]-30)^2)-(xy[1,]-225)*(xy[2,]-30)<=625]_0
w0[((xy[1,]-225)^2+2*(xy[2,]-30)^2)-(xy[1,]-225)*(xy[2,]-30)<=625&xy[2,]>27&xy[2,]<31]_-.5
w0[((xy[1,]-225)^2+2*(xy[2,]-30)^2)-(xy[1,]-225)*(xy[2,]-30)<=625&xy[1,]>223&xy[1,]<227]_.5
w0[((xy[2,]-225)^2+2*(xy[1,]-30)^2)+(xy[2,]-225)*(xy[1,]-30)<=625]_0
w0[((xy[2,]-225)^2+2*(xy[1,]-30)^2)+(xy[2,]-225)*(xy[1,]-30)<=625&xy[1,]>27&xy[1,]<31]_-.5
w0[((xy[2,]-225)^2+2*(xy[1,]-30)^2)+(xy[2,]-225)*(xy[1,]-30)<=625&xy[2,]>223&xy[2,]<227]_.5
w0[((xy[2,]-225)^2+(xy[1,]-225)^2)+1*(xy[2,]-225)*(xy[1,]-225)<=400]_0
w0[((xy[2,]-30)^2+(xy[1,]-30)^2)<=256]_0
sigma_.25
y_w0+rnorm(w0,0,sigma)
yhat_awsbi(y,graph=T,rmax=10)
par(mfrow=c(1,3))
image(y,col=gray((0:255)/255))
image(yhat$yhat,zlim=range(y),col=gray((0:255)/255))
image(w0,zlim=range(y),col=gray((0:255)/255))
}
\keyword{regression}
\keyword{nonparametric}
\keyword{smoothing}

