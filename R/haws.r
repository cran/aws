#
#    R - function  aws  for  Adaptive Weights Smoothing (AWS)
#    in regression models with additive sub-Gaussian errors
#    local constant and local polynomial approach
#
#    Copyright (C) 2002 Weierstrass-Institut f?r
#                       Angewandte Analysis und Stochastik (WIAS)
#
#    Author:  J?rg Polzehl
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
#  USA.
#
awsh <- function(y,x=NULL,p=0,sigma2=NULL,qlambda=NULL,eta=0.5,tau=NULL,
                lkern="Triangle",hinit=NULL,hincr=NULL,hmax=100,hmaxs=2*hmax,
                u=NULL,graph=FALSE,demo=FALSE,symmetric=NULL,conf=FALSE,
                qconf=.95,alpha=2)
{
#
#    first check arguments and initialize
#
args <- match.call()
if(p>0) symmetric <- FALSE
if(is.null(symmetric)) symmetric <- FALSE
if(is.null(qlambda)) {
if(p>5) return("no default for qlambda for p>5")
qlambda <- switch(p+1,.966,.92,.92,.92,.92,.92)
if(symmetric==TRUE) qlambda <- .985
}
if(conf) qconf<-qnorm(1-(1-qconf)/2)
if(qlambda>=1 || qlambda<.6) return("Inappropriate value of qlambda")
if(eta<0 || eta>=1) return("Inappropriate value of eta")
if(demo&& !graph) graph <- TRUE
taudefault <- NULL
# now check which procedure is appropriate
gridded <- is.null(x)
if(gridded){
#  this is the version on a grid
dy <- dim(y)
if(is.null(dy)) {
   form <- "uni"
   ddim  <- 1
   n <- length(y)
   dp1 <- p+1
}
if(length(dy)>1){
return("this is heteroscedastics univariate regression only")
}
} else {
# not gridded
return("currently this is heteroscedastics univariate regression on a grid only")
ddim <- 1
if(is.null(dx)) {
#
#    order data by order of x
#
    form <- "uni"
    n <- length(x)
    if(n!=length(y)) return("incompatible lengths of x and y")
    ox <- order(x)
    x <- x[ox]
    y <- y[ox]
    dp1 <- p+1
}
if(length(y)!=n) return("incompatible dimensions of x and y")
   #
   #
   }
#
#     now set hincr, sigma2 if not provided
#
if(is.null(hincr)) hincr <- 1.25^(1/ddim)
if(is.null(sigma2)){
#
#    heteroscedastic regression estimation of sigma2
#
z<-diff(y)/sqrt(2)
sigma2<-laws(z,hmax=hmaxs,model="Volatility")$theta^2
sigma2<-(sigma2[c(1,1:(n-1))]+sigma2[c(1:(n-1),n-1)])/2
}
if(length(sigma2)==1) sigma2<-rep(sigma2,n)
#
#    now generate kernel on a grid
#
getkern <- function(x,kern)
switch(kern,Triangle=pmax(0,(1-x)),
            Quadratic=pmax(0,(1-x))^2,
            Cubic=pmax(0,(1-x))^3,
            Uniform=as.numeric(abs(x)<=1),
            Exponential=exp(-5*x),
            {
            cat("Triangle kernel is used as default\n");
            pmax(0,(1-x))
            })
# this gives a discretized kern on [0,1.01] for use of (xij^2) as argument
#  length 102  (last element to avoid numerical problems if x(j)==xi+-h)
kernl <- getkern(seq(0,1.01,.01),lkern)
kerns <- getkern(seq(0,1.01,.01),"Exponential")
#
#   get lambda as quantile of appropriate chisq, rescale to be consistent
# with the paper and multiply by 2*sigma2 to get 2*sigma2*lambda in lamakt
#
lamakt <- 10*qchisq(qlambda,dp1)
#
#  set tau in case it may be necessary
#
if(form=="uni") taudefault <- 1.5*3^p else if(is.null(taudefault))
                                          taudefault <- switch(p+1,1,13.5,150)
if(is.null(tau)) tau <- taudefault
#  rescale tau for use in scaled kernel to be consistent with the paper
tau <- 6*tau
#
#    now select the correct aws-procedure
#
#   cases:    gridded      uni   p=0
#             !gridded     uni   p>=0
#             gridded      bi    p=0
#             gridded      bi    p=1,2
#             gridded      tri   p=0
#             !gridded     multi p=0,1
#             !gridded     multi p=0,1  Nearest-Neighbor
#
if(gridded &&  form=="uni" && p==0){
###
###              gridded     uni    p=0
###
###           this should run a little faster than the nongridded version
###
bi <- bi2 <- ai <- theta <- numeric(n)
if(is.null(hinit)||hinit<1) hinit <- 1
#  first initialize
z <- .Fortran("ihawsuni",
              as.double(y),
              as.integer(n),
              as.double(hinit),
              bi=as.double(bi),
              bi2=as.double(bi2),
              ai=as.double(ai),
              as.double(kernl),
              as.double(sigma2),PACKAGE="aws")[c("ai","bi","bi2")]
bi <- z$bi
bi2 <- z$bi2
ai <- z$ai
theta <- ai/bi
if(graph){
par(mfrow=c(1,1),mar=c(3,3,2.5,.5),mgp=c(2,1,0))
plot(y,ylim=range(y,theta),col=3)
if(!is.null(u)) lines(u,col=2)
lines(theta,lwd=2)
if(conf){
lines(theta+qconf/sqrt(bi),col=4)
lines(theta-qconf/sqrt(bi),col=4)
}
title(paste("Reconstruction  h=",signif(hinit,3)))
if(!is.null(u)) cat("bandwidth: ",signif(hinit,3),"   MSE: ",
                    mean((theta-u)^2),"   MAE: ",mean(abs(theta-u)),"\n")
if(demo) readline("Press return")
}
# now run aws-cycle
hakt <- hinit*hincr
if(graph){
#
#   run single steps to display intermediate results
#
while(hakt<=hmax){
z <- .Fortran("lhawsuni",
              as.double(y),
              as.integer(n),
              as.double(hakt),
              as.double(lamakt),
              as.double(theta),
              bi=as.double(bi),
              bi2=as.double(bi),
              ai=as.double(ai),
              as.double(kernl),
              as.double(kerns),
              as.logical(symmetric),
              as.double(sigma2),PACKAGE="aws")[c("ai","bi","bi2")]
ai <- (1-eta)*z$ai + eta * ai
bi <- (1-eta)*z$bi + eta * bi
bi2 <- (1-eta)*z$bi2 + eta * bi2
#  this is correct only if wij = 0 or 1,  or if eta = 0, or if the weights have stabilized
#  but should deliver a reasonable approximation at the end of the iteration process
theta  <- ai / bi
sdtheta<-sqrt(bi2)/bi
plot(y,ylim=range(y,theta),col=3)
if(!is.null(u)) lines(u,col=2)
lines(theta,lwd=2)
if(conf){
lines(theta+qconf*sdtheta,col=4)
lines(theta-qconf*sdtheta,col=4)
}
title(paste("Reconstruction  h=",signif(hakt,3)))
cat("bandwidth: ",signif(hakt,3),
                    "PMSE",mean((theta-y)^2+alpha*bi2/(bi-1/sigma2)^2),mean((theta-y)^2+alpha/(bi-1/sigma2)))
if(!is.null(u)) cat("   MSE: ",
                    mean((theta-u)^2),"   MAE: ",mean(abs(theta-u)))
cat("\n")
if(demo) readline("Press return")
hakt <- hakt*hincr
gc()
}
} else
{
#   run all iterations in one call
z <- .Fortran("ghawsuni",
              as.double(y),
              as.integer(n),
              as.double(hinit),
              as.double(hincr),
              as.double(hmax),
              as.double(lamakt),
              as.double(eta),
              theta=as.double(theta),
              bi=as.double(bi),
              bi2=as.double(bi),
              as.double(ai),
              as.double(kernl),
              as.double(kerns),
              as.double(bi),
              as.logical(symmetric),
              as.double(sigma2),PACKAGE="aws")[c("theta","bi","bi2")]
theta<-z$theta
sdtheta<-sqrt(z$bi2)/z$bi
}
}
      if( form=="uni" && (p>0 || (!gridded && p==0)) ){
###
###                        uni     p>=0
###
if(gridded) x <- 1:length(y)
dp1 <- p+1
dp2 <- p+dp1
bi <- matrix(0,dp2,n)
theta <- ai <- matrix(0,dp1,n)
dxp <- max(diff(x,p+1))*(1+1.e-8)
if(is.null(hinit)||hinit<dxp) hinit <- dxp
#   generate binomial coefficients
cb <- matrix(0,dp1,dp1)
for(i in (1:dp1)) cb[i:dp1,i] <- choose((i:dp1)-1,i-1)
#  first initialize
z <- .Fortran("iphawsun",
              as.integer(n),
              as.integer(dp1),
              as.integer(dp2),
              as.double(x),
              as.double(y),
              as.double(hinit),
              bi=as.double(bi),
              ai=as.double(ai),
              theta=as.double(theta),
              as.double(kernl),
              double(dp1*dp1),
              as.double(sigma2),PACKAGE="aws")[c("ai","bi","theta")]
theta <- matrix(z$theta,dp1,n)
bi <- bi0 <- matrix(z$bi,dp2,n)
ai <- z$ai
if(graph){
par(mfrow=c(1,2),mar=c(3,3,2.5,.5),mgp=c(2,1,0))
plot(x,y,ylim=range(y,theta[1,]),col=3)
if(!is.null(u)) lines(x,u,col=2)
lines(x,theta[1,],lwd=2)
title(paste("Estimated function  h=",signif(hinit,3)))
plot(x,theta[2,],type="l")
title("Estimated first derivative")
}
if(!is.null(u))
cat("bandwidth: ",signif(hinit,3),"   MSE: ",mean((theta[1,]-u)^2),
      "   MAE: ",mean(abs(theta[1,]-u)),"\n")
if(demo) readline("Press return")
# now run aws-cycle
hakt <- hinit*hincr
while(hakt<=hmax){
z <- .Fortran("lphawsun",
              as.integer(n),
              as.integer(dp1),
              as.integer(dp2),
              as.double(x),
              as.double(y),
              as.double(theta),
              as.double(bi),
              bi=as.double(bi),
              as.double(bi0[1,]),
              bi0=as.double(bi0),
              ai=as.double(ai),
              as.double(lamakt),
              as.double(tau),
              as.double(hakt),
              as.double(kernl),
              as.double(kerns),
              as.double(cb),
              double(dp1*dp1),
              double(dp1*dp1),
              double(dp1*dp1),
              double(dp1),
              double(dp1),
              double(dp2),
              double(dp1),
              as.logical(symmetric),
              as.double(sigma2),PACKAGE="aws")[c("ai","bi","bi0")]
    ai <- (1-eta)*z$ai + eta * ai
    bi <- matrix((1-eta)*z$bi + eta * bi,dp2,n)
    bi0 <- (1-eta)*z$bi0 + eta * bi0
    z <- .Fortran("mphawsun",
                  as.integer(n),
                  as.integer(dp1),
                  as.integer(dp2),
                  as.double(ai),
                  as.double(bi),
                  theta=as.double(theta),
                  double(dp1*dp1),
                  sdtheta=double(dp1*n),PACKAGE="aws")[c("theta","sdtheta")]
theta<-matrix(z$theta,dp1,n)
sdtheta<-sqrt(matrix(z$sdtheta,dp1,n))
#  sdtheta contains diaganol elements of the inverse of bi
#  this is correct if weights are 0 or 1, i.e. at the end of the iteration process
if(graph){
plot(x,y,ylim=range(y,theta[1,]),col=3)
if(!is.null(u)) lines(x,u,col=2)
lines(x,theta[1,],lwd=2)
if(conf){
lines(theta[1,]+qconf*sdtheta[1,],col=4)
lines(theta[1,]-qconf*sdtheta[1,],col=4)
}
title(paste("Reconstruction  h=",signif(hakt,3)))
if(conf) ylim<-range(theta[2,]+qconf*sqrt(sdtheta[2,]),theta[2,]-qconf*sqrt(sdtheta[2,])) else ylim<-range(theta[2,])
plot(x,theta[2,],type="l",ylim=ylim)
if(conf){
lines(theta[2,]+qconf*sdtheta[2,],col=4)
lines(theta[2,]-qconf*sdtheta[2,],col=4)
}
title("Estimated first derivative")
}
if(!is.null(u))
cat("bandwidth: ",signif(hakt,3),"   MSE: ",mean((theta[1,]-u)^2),
    "   MAE: ",mean(abs(theta[1,]-u)),"\n")
if(demo) readline("Press return")
hakt <- hakt*hincr
gc()
}
}
###
###            end cases
###
z<-list(theta=theta,sdtheta=sdtheta,y=y,x=x,call=args)
class(z)<-"aws"
z
}
