\name{laws}
\alias{laws}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Likelihood based Adaptive Weights Smoothing }
\description{
  This function implements a Adaptive Weights Smoothing 
  procedure for local constant Poisson, Bernoulli, Exponential,
  Weibull, Volatility and Gaussian models
  as described in Polzehl \& Spokoiny (2002).
}
\usage{
laws(y, x = NULL, qlambda = NULL, eta = 0.5, lkern = "Triangle", 
     model = "Poisson", shape = NULL, hinit = NULL, hincr = NULL, 
     hmax = 10, NN = FALSE, u = NULL, graph = FALSE, demo = FALSE, 
     symmetric = FALSE, wghts=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{y}{ \code{y} contains the observed values at location \code{x}. 
           In case of \code{x==NULL} (second parameter) \code{y} is assumed to be 
           observed on a one, two or three-dimensional grid. The dimension of 
           \code{y} determines if one, two or three-dimensional AWS is used.}
  \item{x}{ \code{x} is either \code{NULL}, in this case \code{y} is assumed 
           to be observed on a grid, or is a matrix, with rows corresponding to variables, 
           containing the design points where \code{y} is observed. }
  \item{qlambda}{ \code{qlambda} determines the scale parameter \code{qlambda} 
           for the stochastic penalty. The scaling parameter in the stochastic 
           penalty \code{lambda} is choosen as the \code{qlambda}-quantile
           of a Chi-square-distribution with number of parameters in the polynomial 
           model as degrees of freedom. If \code{qlambda==NULL} a standard value 
           depending on \code{model} and \code{symmetric} is choosen. }
  \item{eta}{ \code{eta} is a memory parameter used to stabilize the procedure. 
           \code{eta}  has to be between \code{0} and \code{1}, with 
           \code{eta==0.5} being the default. }
  \item{lkern}{ \code{lkern} determines the location kernel to be used. Options 
           are \code{"Uniform"}, \code{"Triangle"}, \code{"Quadratic"}, 
           \code{"Cubic"} and \code{"Exponential"}. Default is \code{"Triangle"}. 
           The Kernel operates on the squared distance, so \code{"Triangle"}
           corresponds to the use of an Epanechnikov kernel in kernel smoothing. 
           \code{"Exponential"} requires larger values of \code{hmax} and 
           therefore more iterations to reach comparable results. }
  \item{model}{ \code{model} determines the distribution type of \code{y}. 
           Currently implemented models are \code{"Poisson"}, \code{"Bernoulli"},
           \code{"Gaussian"}, \code{"Exponential"}, \code{"Weibull"}, \code{"Volatility"}
           (Estimation of the scale parameter of a Gaussian distribution).
           Defaults to \code{"Poisson"}.}
  \item{shape}{ used for additional parameters of the specified distribution if needed,
           i.e. variance if \code{model=="Gaussian"}
            }
  \item{hinit}{ \code{hinit} Initial bandwidth for the location penalty. 
           Appropriate value is choosen in case of \code{hinit==NULL}  }
  \item{hincr}{ \code{hincr} \code{hincr^(1/d)}, with \code{d} the 
           dimensionality of the design, is used as a factor to increase the 
           bandwidth between iterations. Defauts to \code{hincr==1.2}}   
  \item{hmax}{ \code{hmax} Maximal bandwidth to be used. Determines the 
           number of iterations and is used as the stopping rule. }
  \item{NN}{ If \code{NN==TRUE} use nearest neighbor-rules instead of 
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
           \code{sij/lambda}. See references for details.}
   \item{wghts}{ Specifies wghts for distance evaluation on a bi- or trivariate grid.
                Allows for anisotropic local neighborhoods. If \code{wghts=NULL}
                isotropic neighborhoods are used. }
}
\details{ This function implements an adaptive weights 
   smoothing (AWS) procedure for a several classes of distributions for the dependent
   variable in local constant regression models. The approach generalizes the original 
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
    
      
   Currently implemented are the following models (specified by parameter \code{model}):

   \describe{
   \item{Binary response}{\code{model="Bernoulli"}}
   \item{Poisson regression}{\code{model="Poisson"} This model allows e.g. for 
   density estimation or for the analysis of poisson count data on a grid (e.g. Positron
   emission tomography (PET)).}
   \item{Exponential regression}{\code{model="Exponential"} Applications of this model
   include e.g. test for constant (over time) failure rates and estimation of tail indicies.}
   \item{Gaussian regression}{\code{model="Gaussian"} This essentially coincides with the
   local constant regression model with additive sub-gaussian errors as described in 
   Polzehl and Spokoiny (2000, 2003).}
   \item{Weibull regression}{\code{model="Weibull"} Applications in reliability analysis.}
   \item{Volatility model}{\code{model="Volatility"} .}
   }
   
   The essential parameter in the procedure is \code{qlambda}. This parameter has an 
   interpretation as a significance level of a test for equivalence of two local
   parameter estimates. Optimal values mainly depend on the choosen \code{model}
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
\references{ 
\item{ }{Polzehl, J. and Spokoiny, V. (2003). \emph{Varying coefficient 
regression modeling 
     by adaptive weights smoothing}, WIAS-Preprint 818.}
\item{ }{Polzehl, J. and Spokoiny, V. (2002). \emph{Local likelihood 
modelling by adaptive weights smoothing}, WIAS-Preprint 787.}
\item{ }{ Polzehl, J. and Spokoiny, V. (2000). \emph{Adaptive Weights Smoothing
     with applications to image restoration}, J.R.Statist.Soc. B, 62,
     Part 2, pp. 335-354}
 }
\author{ Joerg Polzehl, \email{polzehl@wias-berlin.de}, 
\url{http://www.wias-berlin.de/project-areas/stat/projects/adaptive-image-processing.html}}


\seealso{ SEE ALSO \code{\link{aws}}, \code{\link{awsdens}}, \code{\link{awstindex}} }

\examples{
###
###    Artificial PET data
###
x <- 1:128
f12 <- function(x,y){
2*((1.5*(x-64)^2+(y-64)^2<3025)) +
((x-64)^2+(y-104)^2<16)-((x-92)^2+(y-84)^2<25)+
((x-78)^2+(y-84)^2<25)-((x-64)^2+(y-84)^2<36)+
((x-50)^2+(y-84)^2<36)-((x-36)^2+(y-84)^2<25)+
((x-92)^2+(y-64)^2<25)-((x-78)^2+(y-64)^2<16)+
((x-64)^2+(y-64)^2<16)-((x-50)^2+(y-64)^2<25)+
((x-36)^2+(y-64)^2<25)-((x-92)^2+(y-44)^2<36)+
((x-78)^2+(y-44)^2<36)-((x-64)^2+(y-44)^2<25)+
((x-50)^2+(y-44)^2<25)-((x-36)^2+(y-44)^2<16)+
((x-64)^2+(y-24)^2<16)
}
u0 <- 2*outer(x,x,"f12")
set.seed(1)
y <- matrix(rpois(u0,u0),128,128)
# use larger hmax for good results
yhat <- laws(y,model="Poisson",hmax=4)$theta
par(mfrow=c(1,3),mar=c(3,3,3,.5),mgp=c(2,1,0))
image(y,col=gray((0:255)/255))
title("Observed image")
image(yhat,col=gray((0:255)/255))
title("AWS-Reconstruction")
image(u0,col=gray((0:255)/255))
title("True image")
rm(u0,y,yhat,x)
###
###   VOLATITILTY ESTIMATION
###
###   artificial example
###
sigma <- c(rep(1,125),rep(2,125),rep(.5,125),rep(1,125))
set.seed(1)
x <- rnorm(sigma,0,sigma)
tmpa <- laws(x,model="Volatility",u=sigma,graph=TRUE,hmax=250)
tmps <- laws(x,model="Volatility",u=sigma,hmax=250,symmetric=TRUE)
par(mfrow=c(1,1),mar=c(3,3,3,1))
plot(abs(x),col=3,xlab="time t",ylab=expression(abs( R[t] )))
lines(tmpa$theta,col=1,lwd=2)
lines(tmps$theta,col=1,lwd=2,lty=2)
lines(sigma,col=2,lwd=2,lty=3)
legend(350,5.5,c("asymmetric AWS","symmetric AWS","true sigma"),
                 lwd=c(2,2,2),lty=c(1,2,3),col=c(1,1,2))
title(expression(paste("Estimates of  ",sigma(t))))
rm(tmpa,tmps,sigma,x)
}
\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
\keyword{ regression }% __ONLY ONE__ keyword per line
\keyword{ models }% __ONLY ONE__ keyword per line
