\name{awsuni}
\alias{awsuni}
\title{One-dimensional Adaptive Weights Smoothing}
\description{Performes one-dimensional Adaptive Weigths Smoothing (depreciated version, use aws instead)}
\usage{awsuni(y, lambda=3, gamma=1.3, eta =4, s2hat = NULL, kstar = length(radii),
              radii = c(1:8,(5:12)*2,(7:12)*4,(7:12)*8,(7:10)*16,(6:8)*32,
                      (5:8)*64,(5:8)*128,(5:8)*256),
          rmax=max(radii),graph = FALSE,z0 = NULL, eps = 1e-08, 
          control="dyadic", demomode=FALSE)}
\arguments{
\item{y}{observed values (ordered by value of independent variable)}
\item{lambda}{main smoothing parameter (should be approximately 3)}
\item{gamma}{allow for increase of variances during iteration by factor gamma (!! gamma >=1) }
\item{eta}{main control parameter (should be approximately 4)   }
\item{s2hat}{initial variance estimate (if available,
           can be either a number (homogeneous case), a vector of same length 
           as y (inhomogeneous variance) or NULL (a homogeneous variance estimate
           will be generated in this case)}
\item{kstar}{maximal number of iterations to perform, actual number may be smaller depending
             on parameters radii, rmax and eps}
\item{radii}{radii of neighbourhoods used}
\item{rmax}{maximal radius of neighborhood to be used, may change kstar}
\item{graph}{logical, if TRUE progress (for each iteration) is illustrated grahically,
           if FALSE the program runs until the final estimate is obtained (much faster !!!)}
\item{z0}{allows for submission of "true" values for illustration and test purposes; only
               if graph=TRUE,  MSE and MAE are reported for each iteration step}
\item{eps}{stop iteration if $||(yhatnew - yhat)||^2 < eps * sum(s2hat)$}
\item{control}{the control step is performed in either a dyadic sceme
           ("dyadic") or using all previous estimates (otherwise)}
\item{demomode}{if TRUE the function will wait for user input after each iteration; 
                  only if graph=TRUE }
}

\value{
A list with components
\item{yhat}{estimates of the regression function (corresponding to the y's)}
\item{shat}{estimated standard deviations of yhat (conditional on the chosen weights)}
\item{args}{Main arguments supplied to awsuni}
}
\references{Polzehl, J. and Spokoiny, V. (2000). \emph{Adaptive Weights Smoothing with applications
to image restoration}, J.R.Statist.Soc. B, 62, Part 2, pp. 335-354}
\note{Although the algorithm evaluates a regression model the structure of
the regression function only depends on the ordering of the independent variable.
Therefore no independent variable is to be given as a parameter but the 
values of the dependent variable are required to be ordered by the value
of the independent variable.
This function is superseded by function aws and will be 
removed in the next mayor version of the package.}
\author{Joerg Polzehl \email{polzehl@wias-berlin.de}}
\seealso{\code{\link{aws}}, \code{\link{awsbi}}, \code{\link{awstri}}}

\examples{
#  Blocks data (from Donoho, Johnstone, Kerkyacharian and Picard (1995))
mofx6 <- function(x){
xj <- c(10,13,15,23,25,40,44,65,76,78,81)/100
hj <- c(40,-50,30,-40,50,-42,21,43,-31,21,-42)*.37
Kern <- function(x) (1-sign(x))/2
apply(Kern(outer(xj,x,"-"))*hj,2,sum)
}
x <- seq(0,1,1/2047)
fx6 <- mofx6(x)
#    sigma==3
y <- rnorm(fx6,fx6,3)
tmp <- awsuni(y)
par(mfrow=c(1,1))
plot(x,y)
lines(x,tmp$yhat,col=2)
lines(x,fx6,col=3)
title(expression(paste("AWS Reconstruction of blocks data  ",sigma==3)))
rm(x,y,fx6,mofx6,tmp)
}
\keyword{regression}
\keyword{nonparametric}
\keyword{smooth}
