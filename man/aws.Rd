\name{aws}
\alias{aws}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Local polynomial Adaptive Weights Smoothing for regression with
additive errors  }
\description{ This function implements a local polynomial
              Adaptive Weights Smoothing procedure for
              regression problems with additive errors as
              described in Polzehl \& Spokoiny (2003).
              This function superseedes functions awsuni, awsbi and awstri.
}
\usage{
aws(y, x = NULL, p = 0, sigma2 = NULL, qlambda = NULL, eta = 0.5, 
    tau = NULL, lkern = "Triangle", hinit = NULL, hincr = NULL, 
    hmax = 10, NN = FALSE, u = NULL, graph = FALSE, demo = FALSE, 
    symmetric = NULL, wghts= NULL)
}
\arguments{
  \item{y}{ \code{y} contains the observed values (regression function plus errors).
           In case of \code{x=NULL} (second parameter) \code{y} is assumed to be
           observed on a one, two or three-dimensional grid. The dimension of
           \code{y} determines if one, two or three-dimensional AWS is used.}
  \item{x}{ \code{x} is either \code{NULL}, in this case \code{y} is assumed
           to be observed on a grid, or is a matrix, with rows corresponding to variables,
           containing the design points where \code{y} is observed. }
  \item{p}{ \code{p} is the degree of the polynomial model to use. For univariate
           regression \code{p} can be an nonnegative integer less or equal than 5, for
           bivariate regression on a grid \code{p} can be either \code{0}, \code{1} or
           \code{2}, for three dimensional AWS only \code{p=0} is implemented.  }
  \item{sigma2}{ \code{sigma2} can be used to provide an estimate for the error
           variance. If \code{is.null(sigma2)} a variance estimate is generated from
           the data.}
  \item{qlambda}{ \code{qlambda} determines the scale parameter \code{qlambda}
           for the stochastic penalty. The scaling parameter in the stochastic
           penalty \code{lambda} is choosen as the \code{qlambda}-quantile
           of a Chi-square-distribution with number of parameters in the polynomial
           model as degrees of freedom. If \code{is.null(qlambda)} a standard value
           depending on \code{p} is choosen.}
  \item{eta}{ \code{eta} is a memory parameter used to stabilize the procedure.
           \code{eta}  has to be between \code{0} and \code{1}, with
           \code{eta=.5} being the default. }
  \item{tau}{ \code{tau} is used in case of a the scale parameter polynomial
           degree \code{p!=0} only. It is the scale parameter in the extention
           penalty used to prevent from leverage problems. The default value
           for \code{tau} depends on \code{p} and on the number of parameters
           in the local model. }
  \item{lkern}{ \code{lkern} determines the location kernel to be used. Options
           are \code{"Uniform"}, \code{"Triangle"}, \code{"Quadratic"},
           \code{"Cubic"} and \code{"Exponential"}. Default is \code{"Triangle"}.
           The Kernel operates on the squared distance, so \code{"Triangle"}
           corresponds to the use of an Epanechnikov kernel in kernel smoothing.
           \code{"Exponential"} requires larger values of \code{hmax} and
           therefore more iterations to reach comparable results. }
  \item{hinit}{ \code{hinit} Initial bandwidth for the location penalty.
           Appropriate value is choosen in case of \code{hinit=NULL} }
  \item{hincr}{ \code{hincr} \code{hincr^(1/d)}, with \code{d} the
           dimensionality of the design, is used as a factor to increase the 
           bandwidth between iterations. Defauts to \code{hincr=1.25} }
  \item{hmax}{ \code{hmax} Maximal bandwidth to be used. Determines the
           number of iterations and is used as the stopping rule. }
   \item{NN}{ If \code{NN=TRUE} use nearest neighbor-rules instead of
           distances in the location term. }
  \item{u}{ \code{u} used to supply values of the true regression function 
           for test purposes to calculate  Mean Squared Error (MSE) and
           Mean Absolute Error (MAE) }
  \item{graph}{ \code{graph} if \code{TRUE} results are displayed after each
           iteration step. }
  \item{demo}{ \code{demo} if \code{TRUE} after each iteration step results
           are displayed and the process waits for user interaction. }
  \item{symmetric}{ If \code{symmetric==TRUE} the stochastic penalty is
           symmetrized, i.e. \code{(sij + sji)/lambda} is used instead of
           \code{sij/lambda}. See references for details.
           \code{symmetric==FALSE} is forced if \code{p!=0}}
   \item{wghts}{ Specifies wghts for distance evaluation on a bi- or trivariate grid.
                Allows for anisotropic local neighborhoods. If \code{wghts=NULL}
                isotropic neighborhoods are used. }
}
\details{This function implements an adaptive weights smoothing (AWS) procedure for
   local polynomial models with additive errors. The approach generalizes the original
   AWS procedure from Polzehl and Spokoiny (2000).

   Adaptive weights smoothing is an iterative data adaptive smoothing technique that
   is designed  for smoothing in regression problems with discontinuous regression
   function. The basic assumption is that the regression function can be approximated
   by a simple local,  e.g. local constant or local polynomial, model.
   The estimate of the regression function, i.e. the conditional expectation of \code{y}
   given \code{x} is computed as a weighted maximum likelihood estimate, with weights choosen
   in a completely data adaptive way. The procedure is edge preserving. If the assumed local
   model is globally valid, almost all weights used will be 1, i.e. the resulting estimate
   almost is the global estimate.


   Currently implemented are the following models (specified by parameter \code{p} and attributes
   of \code{x} and \code{y}) are implemented:

   \describe{
   \item{p=0, x=NULL}{local constant univariate, bivariate or three-variate AWS on a
             grid specified by the dimensionality of \code{y}. Measurement unit for \code{hinit}
             and \code{hmax} is 1 (number grid points within distance \code{h})}
   \item{p=0, !is.null(x)}{local constant AWS for arbitrary design given in \code{x}. If
             \code{x} is a matrix, then the rows correspond to components (variables).
             If \code{NN==TRUE} distances in the location penalty are replaced by nearest-
             neighbor-rules. }
   \item{p=1 or 2, x=NULL}{local linear univariate and bivariate AWS on a grid
             grid specified by the dimensionality of \code{y}. Measurement unit for \code{hinit}
             and \code{hmax} is 1 (number grid points within distance \code{h})}
   \item{p=1, !is.null(x)}{local linear AWS for arbitrary design given in \code{x}. If
             \code{x} is a matrix, then the rows correspond to components (variables).
             If \code{NN==TRUE} distances in the location penalty are replaced by nearest-
             neighbor-rules. }
   \item{\code{p>2}, univariate}{univariate local polynomial AWS for arbitrary univariate design
              given in \code{x}. \code{y} has to be a vector. If
             \code{x} is provided it has to be a vector of the same length as \code{y}.
             If \code{is.null(x)} an equidistant design (grid) with grid-step 1 is assumed.}
         }

   The essential parameter in the procedure is \code{qlambda}. This parameter has an
   interpretation as a significance level of a test for equivalence of two local
   parameter estimates. Optimal values mainly depend on the choosen \code{p}
   and the value of \code{symmetric} (determines the use of an asymmetric or a symmetrized
   test). The optimal values only slightly depend on the model parameters, i.e. the
   default parameters should work in most situations. Larger values of \code{qlambda}
   may lead to oversmoothing, small values of \code{qlambda} lead to a random segmentation
   of homogeneous regions. A good value of \code{qlambda} can be obtained by the propagation
   condition, requiring that in case of global validity of the local model the 
   estimate for large \code{hmax} should be equal to the global estimate.

   The numerical complexity of the procedure is mainly determined by \code{hmax}. The number
   of iterations is \code{d*log(hmax/hinit)/log(hincr)} with \code{d} being the dimension
   of \code{y}. Comlexity in each iteration step is \code{Const*hakt*n} with \code{hakt}
   being the actual bandwith in the iteration step and \code{n} the number of design points.
   \code{hmax} determines the maximal possible variance reduction.
}
\value{
  \item{theta }{Parameter estimates, first dimension corresponds to parameter
components}
  \item{y}{values provided in \code{y}}
  \item{x}{values provided in \code{x}}
  \item{call}{actual function call}
}
\references{ \item{ }{ Polzehl, J. and Spokoiny, V. (2003).
   \emph{Varying coefficient regression modeling
     by adaptive weights smoothing}, WIAS-Preprint 818. }
\item{ }{ Polzehl, J. and Spokoiny, V. (2000). \emph{Adaptive Weights Smoothing
     with applications to image restoration}, J.R.Statist.Soc. B, 62,
     Part 2, pp.335-354
 }
 }
\author{ Joerg Polzehl, \email{polzehl@wias-berlin.de},
\url{http://www.wias-berlin.de/project-areas/stat/projects/adaptive-image-processing.html}
 }
\note{ This function supersedes functions \code{\link{awsuni}}, \code{\link{awsbi}}, and \code{\link{awstri}} }


\seealso{ SEE ALSO \code{\link{laws}}, \code{\link{awsuni}}, \code{\link{awsbi}}, 
\code{\link{awstri}} }

\examples{
#######################################################
#                                                     
#                  univariate examples                
#                                                     
####################################################### 
#  
#    Blocks data (Example 6 from Fan & Gijbels (1996)
#
     mofx6 <- function(x){
     xj <- c(10,13,15,23,25,40,44,65,76,78,81)/100
     hj <- c(40,-50,30,-40,50,-42,21,43,-31,21,-42)*.37
     Kern <- function(x) (1-sign(x))/2
     apply(Kern(outer(xj,x,"-"))*hj,2,sum)
     }
#
#    sigma==1   step by step with graphics
#
     fx6 <- mofx6(seq(0,1,1/2047))
     y <- rnorm(fx6,fx6,1)
     tmp <- aws(y,p=0,hmax=100,graph=TRUE)
     par(mfrow=c(1,1),mar=c(3,3,2.5,.5),mgp=c(2,1,0))
     plot(seq(0,1,1/2047),y)
     lines(seq(0,1,1/2047),tmp$theta,col=3,lwd=2)
     lines(seq(0,1,1/2047),fx6,col=2)
#
#    sigma==3   without graphics (much faster)
#
     y <- rnorm(fx6,fx6,3)
     tmp <- aws(y,hmax=500)
     par(mfrow=c(1,1),mar=c(3,3,2.5,.5),mgp=c(2,1,0))
     plot(seq(0,1,1/2047),y)
     lines(seq(0,1,1/2047),tmp$theta,col=3,lwd=2)
     lines(seq(0,1,1/2047),fx6,col=2)
     rm(mofx6,fx6,y,tmp)
#  
#    second example from Polzehl and Spokoiny (2002)
#
     f2 <- function(x) sin(2*pi*1.2/(x+.2))
     n <- 1000
     x <- seq(0,1,length=n)
     fx2 <- f2(x)
     set.seed(1)
     sigma <- .25
     y <- rnorm(x,fx2,sigma)
#   increase hmax to 2 for good results
     ex1p0 <- aws(y,x,p=0,hmax=.1)$theta[1,]
     ex1p1 <- aws(y,x,p=1,hmax=.1)$theta[1,]
     ex1p2 <- aws(y,x,p=2,hmax=.1)$theta[1,]
     ex1p3 <- aws(y,x,p=3,hmax=.1)$theta[1,]
     par(mfrow=c(2,2),mar=c(2.5,2.5,2.5,.5),mgp=c(1.5,.5,0))
     plot(x,y)
     lines(x,fx2,col=2)
     lines(x,ex1p0,col=3,lwd=2)
     title("local constant AWS")
     plot(x,y)
     lines(x,fx2,col=2)
     lines(x,ex1p1,col=3,lwd=2)
     title("local linear AWS")
     plot(x,y)
     lines(x,fx2,col=2)
     lines(x,ex1p2,col=3,lwd=2)
     title("local quadratic AWS")
     plot(x,y)
     lines(x,fx2,col=2)
     title("local cubic AWS")
     lines(x,ex1p3,col=3,lwd=2)
     rm(f2,fx2,x,sigma,y,ex1p0,ex1p1,ex1p2,ex1p3)
#######################################################
#                                                     
#                  bivariate examples                 
#                                                     
####################################################### 
#  
#  Local constant image from Polzehl and Spokoiny (2000)
#
     xy <- rbind(rep(0:255,256),rep(0:255,rep(256,256)))
     indw <- c(1:12,29:48,73:100,133:168,209:256)
     w0 <- matrix(rep(FALSE,256*256),ncol=256)
     w0[indw,] <- TRUE
     w0[,indw] <- !w0[,indw]
     w0 <- w0-.5
     
     w0[((xy[1,]-129)^2+(xy[2,]-129)^2)<=10000&((xy[1,]-129)^2+(xy[2,]-129)^2)>=4900] <- 0
     w0[abs(xy[1,]-xy[2,])<=20&((xy[1,]-129)^2+(xy[2,]-129)^2)<4900] <- 0
     w0[((xy[1,]-225)^2+2*(xy[2,]-30)^2)-(xy[1,]-225)*(xy[2,]-30)<=625] <- 0
     
     w0[((xy[1,]-225)^2+2*(xy[2,]-30)^2)-(xy[1,]-225)*(xy[2,]-30)<=625&xy[2,]>27&xy[2,]<31] <- -.5
     
     w0[((xy[1,]-225)^2+2*(xy[2,]-30)^2)-(xy[1,]-225)*(xy[2,]-30)<=625&xy[1,]>223&xy[1,]<227] <- .5
     w0[((xy[2,]-225)^2+2*(xy[1,]-30)^2)+(xy[2,]-225)*(xy[1,]-30)<=625] <- 0
     
     w0[((xy[2,]-225)^2+2*(xy[1,]-30)^2)+(xy[2,]-225)*(xy[1,]-30)<=625&xy[1,]>27&xy[1,]<31] <- -.5
     
     w0[((xy[2,]-225)^2+2*(xy[1,]-30)^2)+(xy[2,]-225)*(xy[1,]-30)<=625&xy[2,]>223&xy[2,]<227] <- .5
     w0[((xy[2,]-225)^2+(xy[1,]-225)^2)+1*(xy[2,]-225)*(xy[1,]-225)<=400] <- 0
     w0[((xy[2,]-30)^2+(xy[1,]-30)^2)<=256] <- 0
     rm(xy,indw)
     set.seed(1)
     sigma <- .5
     y <- w0+rnorm(w0,0,sigma)
#    run some steps interactively  (increase hmax)
     tmp <- aws(y,graph=TRUE,hmax=2,demo=TRUE)
#    now without graphics and larger hmax
#   increase hmax for better results
     tmp <- aws(y,hmax=5)
     par(mfrow=c(1,3))
     image(y,col=gray((0:255)/255))
     image(tmp$theta,zlim=range(y),col=gray((0:255)/255))
     image(w0,zlim=range(y),col=gray((0:255)/255))
     rm(y,w0,tmp,sigma)
#
#   piecewise smooth image from Polzehl and Spokoiny (2003)
#
     x <- (0:99)/99
     fbi <- function(x,y) (x^2+y^3)*sign(x^2-1*x*y-.75*y^3)
     z0 <- outer(2*(x-.5),2*(x-.25),FUN=fbi)
     z <- z0+rnorm(z0,0,.5)
     par(mfrow=c(1,3),mar=c(3,3,3,.5),mgp=c(2,1,0))
     set.seed(1)
     persp(x,x,z0,phi=15,theta=150,d=5,r=2,expand=1.5,col="white")
     title("True function")
     persp(x,x,z,phi=15,theta=150,d=5,r=2,expand=1.5,col="white")
     title("Observed values")
     image(x,x,z,col=gray((0:255)/255))
     title("Observed values")
#
#   local constant
#
     ex1p0 <- aws(z,hmax=3,symmetric=FALSE)$theta # use hmax=5 or larger
#
#   local linear
#
     ex1p1 <- aws(z,p=1,hmax=3)$theta[1,,]  # use hmax=12 or larger
#
#   local quadratic
#
     ex1p2 <- aws(z,p=2,hmax=3)$theta[1,,]  # use hmax=20 or larger
#
#   display results
#
     par(mfrow=c(2,2),mar=c(0.25,.5,.5,.5),mgp=c(.5,.5,0))
     persp(x,x,z,phi=15,theta=150,d=5,r=2,expand=1,xlab="",ylab="",zlab="",box=FALSE)
     title("original data",line=-3)
     persp(x,x,ex1p0,phi=15,theta=150,d=5,r=2,expand=1,xlab="",ylab="",zlab="",box=FALSE)
     title("local constant AWS",line=-3)
     persp(x,x,ex1p1,phi=15,theta=150,d=5,r=2,expand=1,xlab="",ylab="",zlab="",box=FALSE)
     title("local linear AWS (small hmax)",line=-3)
     persp(x,x,ex1p2,phi=15,theta=150,d=5,r=2,expand=1,xlab="",ylab="",zlab="",box=FALSE)
     title("local quadratic AWS (small hmax)",line=-3)
     rm(x,fbi,z0,ex1p0,ex1p1,ex1p2)

#######################################################
#
#                  three-variate examples
#
#######################################################
#
#  Local constant image from Polzehl and Spokoiny (2000)
#
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
     rm(xy)
     sigma <- .5
     set.seed(1)
     y <- w3+rnorm(w3,0,sigma)
#   increase hmax for reasonable results (hmax >=5)
     tmp <- aws(y,hmax=2)
     par(mfrow=c(1,3))
     for(i in 14:18){
     image(y[,,i],col=gray((0:255)/255))
     image(tmp$theta[,,i],zlim=range(y),col=gray((0:255)/255))
     image(w3[,,i],zlim=range(y),col=gray((0:255)/255))
#     readline()
     }
     rm(w3,y,tmp,sigma)
}
\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
\keyword{ regression }% __ONLY ONE__ keyword per line
\keyword{ models }% __ONLY ONE__ keyword per line
