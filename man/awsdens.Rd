\name{awsdens}
\alias{awsdens}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ local constant density estimation }
\description{
  This function uses adaptive weights smoothing for density estimation. The estimate
  is constructed using poisson regression and the structural assumption of a local
  constant density. Implementation is for 1D, 2D and 3D problems, although it 
  applications may be interesting only in 2D and 3D.
}
\usage{
awsdens(y, ngrid = NULL, nempty = NULL, qlambda = NULL, eta = 0.5, 
        lkern = "Triangle", fu = NULL, hinit = 1, hincr = 1.2, 
    hmax = NULL, graph = FALSE, demo = FALSE, symmetric = TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{y}{ \code{y} contains the observed random variables or vectors. 
       In case of random vectors \code{y} is supposed to be an
       observation matrix with columns containing the observed vectors}
  \item{ngrid}{ \code{ngrid} determines the size of the grid used for binning. 
       If \code{ngrid} is a vector its components are supposed to determine 
       the number of intervals for the corresponding component of \code{y};
       if it is a scalar \code{rep(ngrid,dim(y)[1])} is used.
       If \code{is.null(ngrid)} a default of \code{2*n^{1/d}}, with 
       \code{dim(y)=c(d,n)}, is chosen, i.e. a grid with \code{2^d*n} bins.}
  \item{nempty}{ determines the width of a boundary region of the grid containing only empty bins.
        \code{nempty} defaults to \code{0.1*ngrid}. The grid generated is equispaced with 
        \code{prod(ngrid)} bins and grid-coordinates \code{(nempty+1)[i]} and \code{(ngrid+1-nempty)[i]}
        determined by \code{min(y[i,])} and \code{max(y[i,])}.}
  \item{qlambda}{ \code{qlambda} determines the scale parameter \code{qlambda} 
           for the stochastic penalty. The scaling parameter in the stochastic 
           penalty \code{lambda} is choosen as the \code{qlambda}-quantile
           of a Chi-square-distribution with number of parameters in the polynomial 
           model as degrees of freedom. If \code{qlambda=NULL} a standard value 
           depending on \code{model} and \code{symmetric} is choosen. }
  \item{eta}{ \code{eta} is a memory parameter used to stabilize the procedure. 
           \code{eta}  has to be between \code{0} and \code{1}, with 
           \code{eta=.5} being the default. }
  \item{lkern}{ \code{lkern} determines the location kernel to be used. Options 
           are \code{"Uniform"}, \code{"Triangle"}, \code{"Quadratic"}, 
           \code{"Cubic"} and \code{"Exponential"}. Default is \code{"Triangle"}. 
           The Kernel operates on the squared distance, so \code{"Triangle"}
           corresponds to the use of an Epanechnikov kernel in kernel smoothing. 
           \code{"Exponential"} requires larger values of \code{hmax} and 
           therefore more iterations to reach comparable results. }
  \item{fu}{ \code{fu} can be used to specify a function to calculate the 
           values of the true density on the grid. This may be used
           for test purposes to calculate  Mean Squared Error (MSE) and 
           Mean Absolute Error (MAE) }
  \item{hinit}{ \code{hinit} Initial bandwidth for the location penalty. 
           Appropriate value is choosen in case of \code{hinit==NULL}  }
  \item{hincr}{ \code{hincr} \code{hincr^(1/d)}, with \code{d} the 
           dimensionality of the design, is used as a factor to increase the 
           bandwidth between iterations. Defauts to \code{hincr=1.2} }  
  \item{hmax}{ \code{hmax} Maximal bandwidth to be used. Determines the 
           number of iterations and is used as the stopping rule. }
  \item{graph}{ \code{graph} if \code{TRUE} results are displayed after each 
           iteration step. }
  \item{demo}{ \code{demo} if \code{TRUE} after each iteration step results 
           are displayed and the process waits for user interaction. }
  \item{symmetric}{ If \code{symmetric==TRUE} the stochastic penalty is
           symmetrized, i.e. \code{(sij + sji)/lambda} is used instead of 
           \code{sij/lambda}. See references for details.}
}
\details{
  The density estimate ist constructed using an approach proposed by 
  Lindsay (1974a,b). 1, 2 or 3-dimensional binning is used to produce
  counts, i.e. numbers of observations, per bin, with bins defined by
  a regular grid. THe bin-counts are then considered as poisson random
  variables with intensity depending on the location within the grid.
  Adaptive Weights Smoothing for poisson regression models, i.e. function 
  \code{\link{laws}} with parameter \code{model="Poisson"}, is used 
  to construct a location dependent intensity estimate which is then 
  transformed into a density estimate, see Section 6 in Polzehl and Spokoiny
  (2002b) for details.
}
\value{ The returned object is a list with components
  \item{bin }{Bin counts}
  \item{dens }{density values for the bin. Values correspond to the grid-centers
              determined by component \code{xgrid}}
  \item{xgrid }{list with components containing the coordinates of bin-centers.
            The \code{dim(y)[1]} components of \code{xgrid} correspond to the
            components of the grid-coordinates. }
  \item{call }{actual function call}
}
\references{ 
\item{ }{    Polzehl, J. and Spokoiny, V. (2002). 
\emph{Local likelihood modelling by adaptive weights smoothing}, 
WIAS-Preprint 787 } 
   \item{ }{  Polzehl, J. and Spokoiny, V. (2000). \emph{Adaptive Weights Smoothing
     with applications to image restoration}, J.R.Statist.Soc. B, 62,
     Part 2, pp.335-354}
\item{ }{
Lindsay, J. (1974a). \emph{Comparison of probabbility distributions}, J. Royal Statist. Soc. Ser. B
36, 38-47.}
\item{ }{
Lindsay, J. (1974b). \emph{Construction and comparison of statistical models}, J. Royal Statist. Soc. Ser. B
36, 418-425.
}
}
\author{ Joerg Polzehl, \email{polzehl@wias-berlin.de }}


\seealso{ SEE ALSO \code{\link{aws}}, \code{\link{laws}} }

\examples{
###
###   univariate density estimation
###
f0 <- function(x) 1.5*(x>=0)-(x>=.25)+(x>=.75)-1.5*(x>1)
set.seed(1)
m <- 1000
x1 <- runif(m)
ind <- (runif(m)<f0(x1)/1.5)
n <- 200
y <- x1[ind][1:n]
tmp0 <- awsdens(y,440,20,hmax=2000)
plot(tmp0$xgrid[[1]],tmp0$dens,type="l",lwd=2,
     ylim =range(0,1.5,tmp0$dens),xlab="x",ylab="Density")
lines(tmp0$xgrid[[1]],f0(tmp0$xgrid[[1]]),lty=3,col=2,lwd=2)
title("Estimated and true density (n=200)")
###
###   bivariate density estimation
###
f1 <- function(x,y) 7.5*pmax(x*(1-x^2-y^2),0)*(x>0)*(y>0)
set.seed(1)
m <- 10000
x1 <- runif(m)
x2 <- runif(m)
fx12 <- f1(x1,x2)
ind <- (runif(m)<f1(x1,x2)/.385/7.5)
n <- 2500
y <- rbind(x1[ind],x2[ind])[,1:n]
tmp <- awsdens(y,120,10,hmax=5)
image(tmp$xgrid[[1]],tmp$xgrid[[2]],tmp$dens,xlab="x1",ylab="x2",col=gray(.5+(255:0)/511),lwd=2)
#   thats the estimated density on the grid
lines(seq(0,1,.01),sqrt(1-seq(0,1,.01)^2),col=2,lty=2,lwd=2)
lines(c(0,1),c(0,0),col=2,lty=2,lwd=2)
lines(c(0,0),c(0,1),col=2,lty=2,lwd=2)
#    thats the boundary of the support
title("Estimated density (n=2500)")
#    now add contour lines of the estimated density
contour(tmp$xgrid[[1]],tmp$xgrid[[2]],tmp$dens,nlevels=20,add=TRUE,col=1,lty=1,labcex=.1)
rm(f0,m,x1,ind,n,y,tmp0,f1,x2,fx12,tmp)
}
\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
\keyword{ models }% __ONLY ONE__ keyword per line
