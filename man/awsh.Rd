\name{awsh}
\alias{awsh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Univariate local polynomial Adaptive Weights Smoothing for regression with
heteroscedastic additive errors }
\description{
  This function implements a univariate local polynomial
              Adaptive Weights Smoothing procedure for
              regression problems with heteroskedastic additive errors as
              described in Polzehl \& Spokoiny (2003).
}
\usage{
awsh(y, x = NULL, p = 0, sigma2 = NULL, qlambda = NULL, eta = 0.5, tau = NULL,
lkern = "Triangle", hinit = NULL, hincr = NULL, hmax = 100, hmaxs= 2*hmax, u = NULL,
graph = FALSE, demo = FALSE, symmetric = NULL, conf = FALSE, qconf = 0.95, alpha = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ \code{y} contains the observed values (regression function plus errors).
           In case of \code{x=NULL} (second parameter) \code{y} is assumed to be
           observed on a one-dimensional grid. }
  \item{x}{ \code{x} is either \code{NULL}, in this case \code{y} is assumed
           to be observed on a grid, or is a vector determining the design. }
  \item{p}{ \code{p} is the degree of the polynomial model to use. For univariate
           regression \code{p} can be an nonnegative integer less or equal than 5.  }
  \item{sigma2}{ \code{sigma2} can be used to provide an estimate for the error
           variance. If \code{is.null(sigma2)} a homoskedastic model is assumed and
           a variance estimate is generated from the data. If \code{length(sigma2)==length(y)}
	   this parameter provides variances at the design points.}
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
  \item{hincr}{ \code{hincr} \code{hincr} is used as a factor to increase the
           bandwidth between iterations. Defauts to \code{hincr=1.25} }
  \item{hmax}{ \code{hmax} Maximal bandwidth to be used. Determines the
           number of iterations and is used as the stopping rule. }
  \item{hmaxs}{ \code{hmaxs} Maximal bandwidth to be used when estimating the
           heterogenous variance from consequtive differences of \code{y} by function
           \code{laws}. Determines the number of iterations of \code{laws}. }
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
  \item{conf}{ \code{conf} if \code{TRUE} conditional (on weights) confidence intervals are provided
           at each design point.}
  \item{qconf}{ \code{qconf} determines the level of the conditional (on weights) confidence intervals }
  \item{alpha}{ \code{alpha} Parameter used for a penalized MSE estimate for \code{p=0}. This is
                experimental to try to select \code{hmax}. }
}
\details{
   This function implements an local polynomial adaptive weights smoothing (AWS) procedure for a
   univariate heteroskedastic regression model. The approach generalizes the original
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
   \item{p=0, x=NULL}{local constant univariate AWS on a grid. Measurement unit for \code{hinit}
             and \code{hmax} is 1 (number grid points within distance \code{h})}
   \item{p=0, !is.null(x)}{local constant AWS for arbitrary design given in \code{x}. }
   \item{\code{p>2}}{univariate local polynomial AWS for arbitrary univariate design
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
  \item{theta }{Parameter estimates, first dimension corresponds to parameter components}
  \item{theta }{Standard deviations (conditional on weights) of parameter estimates,
  first dimension corresponds to parameter components}
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
\seealso{  \code{\link{aws}} }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
\keyword{ regression }% __ONLY ONE__ keyword per line
\keyword{ models }% __ONLY ONE__ keyword per line
