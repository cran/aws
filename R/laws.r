#
#    R - function  aws  for  Adaptive Weights Smoothing (AWS)
#    in regression models with additive sub-Gaussian errors               
#    local constant and local polynomial approach                         
#
#    Copyright (C) 2002 Weierstrass-Institut für                          
#                       Angewandte Analysis und Stochastik (WIAS)         
#
#    Author:  Jörg Polzehl                                                
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
aws <- function(y,x=NULL,p=0,sigma2=NULL,qlambda=NULL,eta=0.5,tau=NULL,
                lkern="Triangle",hinit=NULL,hincr=NULL,hmax=10,NN=FALSE,
                u=NULL,graph=FALSE,demo=FALSE,symmetric=NULL,wghts=NULL)
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
if(qlambda>=1 || qlambda<.6) return("Inappropriate value of qlambda")
if(eta<0 || eta>=1) return("Inappropriate value of eta")
if(demo&& !graph) graph <- TRUE
taudefault <- NULL
# now check which procedure is appropriate
gridded <- is.null(x)
if(gridded){
#  this is the version on a grid
if(is.null(hinit)||hinit<=0) hinit <- 1
dy <- dim(y)
if(is.null(dy)) {
   form <- "uni"
   ddim  <- 1
   n <- length(y)
   dp1 <- p+1
}
if(length(dy)==2){
   form <- "bi"
   ddim  <- 2
if(is.null(wghts)) wghts<-c(1,1)
hinit<-hinit/wghts[1]
hmax<-hmax/wghts[1]
wghts<-(wghts[2]/wghts[1])^2
#  only use a wght for the second component
n1 <- dy[1]
n2 <- dy[2]
n <- n1*n2
if(p>2) return("bivariate aws on a grid is not implemented for p>2")
dp1 <- switch(p+1,1,3,6)
}
if(length(dy)==3){
   form <- "tri"
   ddim  <- 3
if(is.null(wghts)) wghts<-c(1,1,1)
hinit<-hinit/wghts[1]
hmax<-hmax/wghts[1]
wghts<-(wghts[2:3]/wghts[1])^2
#  only use a wght for the second and third component
n1 <- dy[1]
n2 <- dy[2]
n3 <- dy[3]
n <- n1*n2*n3
dp1 <- 3*p+1
}
if(length(dy)>3)
   return("AWS for more than 3 dimensional grids is not implemented")
} else {
# not gridded
dx <- dim(x)
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
}else {
   px <- dx[1]
   n <- dx[2]
   if(p>1) {
      p <- 1
      cat("p is set to 1, the maximal polynomial degree implemented")
      }
   form <- "multi"
   if(n!=length(y)) return("incompatible dimensions of x and y")
   if(is.null(wghts)||length(wghts)!=px) wghts <- rep(1,px)
   dp1 <- 1+p*px
   taudefault <- 1.5*3^px
#
#  now generate matrix of nearest neighbors
#  hmax is interpreted as maximal number of neighbors
#
   if(NN){
   ihmax <- trunc(hmax)
   if(ihmax>n) ihmax <- n
   neighbors <- matrix(0,ihmax,n)
   for (i in 1:n) {
      if(px==1) adist <- (x-x[i])^2 else adist <- wghts%*%((x-x[,i])^2)
      neighbors[,i] <- order(adist)[1:ihmax]
      }
   } else {
   ihmax <- n
   ddim <- px
   neighbors <- distmat <- matrix(0,n,n)
   for (i in 1:n) {
      if(px==1) adist <- (x-x[i])^2 else adist <- wghts%*%((x-x[,i])^2)
      od <- order(adist)
      distmat[,i] <- adist[od]
      neighbors[,i] <- od
      }
#  now reduce memory used to whats needed                                 
   gc()
   distmat <- sqrt(distmat)
   maxdist <- apply(distmat,1,max)
   ihmax <- sum(maxdist<=hmax)
   distmat <- distmat[1:ihmax,]
   neighbors <- neighbors[1:ihmax,]
   maxdist <- maxdist[1:ihmax]
   gc()
   }
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
IQRdiff <- function(y) IQR(diff(y))/1.908
sigma2 <- IQRdiff(y)^2
}
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
lamakt <- 10*qchisq(qlambda,dp1)*sigma2
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
bi <- ai <- theta <- numeric(n)
#  first initialize
z <- .Fortran("iawsuni",
              as.double(y),
              as.integer(n),
              as.double(hinit),
              bi=as.double(bi),
              ai=as.double(ai),
              as.double(kernl),PACKAGE="aws")[c("ai","bi")]
bi <- z$bi
ai <- z$ai
theta <- ai/bi
if(graph){
par(mfrow=c(1,2),mar=c(3,3,2.5,.5),mgp=c(2,1,0))
plot(y,ylim=range(y,theta),col=3)
if(!is.null(u)) lines(u,col=2)
lines(theta,lwd=2)
title(paste("Reconstruction  h=",signif(hinit,3)))
plot(bi,type="l")
title("Sum of weights")
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
z <- .Fortran("lawsuni",
              as.double(y),
              as.integer(n),
              as.double(hakt),
              as.double(lamakt),
              as.double(theta),
              bi=as.double(bi),
              ai=as.double(ai),
              as.double(kernl),
              as.double(kerns),
              as.logical(symmetric),PACKAGE="aws")[c("ai","bi")]
ai <- (1-eta)*z$ai + eta * ai
bi <- (1-eta)*z$bi + eta * bi
theta  <- ai / bi
plot(y,ylim=range(y,theta),col=3)
if(!is.null(u)) lines(u,col=2)
lines(theta,lwd=2)
title(paste("Reconstruction  h=",signif(hakt,3)))
plot(bi,type="l")
title("Sum of weights")
if(!is.null(u)) cat("bandwidth: ",signif(hakt,3),"   MSE: ",
                    mean((theta-u)^2),"   MAE: ",mean(abs(theta-u)),"\n")
if(demo) readline("Press return")
hakt <- hakt*hincr
gc()
}
} else
{
#   run all iterations in one call
theta <- .Fortran("gawsuni",
              as.double(y),
              as.integer(n),
              as.double(hinit),
              as.double(hincr),
              as.double(hmax),
              as.double(lamakt),
              as.double(eta),
              theta=as.double(theta),
              as.double(bi),
              as.double(ai),
              as.double(kernl),
              as.double(kerns),
              as.double(bi),
              as.logical(symmetric),PACKAGE="aws")$theta
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
z <- .Fortran("ipawsuni",
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
              double(dp1*dp1),PACKAGE="aws")[c("ai","bi","theta")]
theta <- matrix(z$theta,dp1,n)
bi <- bi0 <- matrix(z$bi,dp2,n)
ai <- z$ai
if(graph){
par(mfrow=c(1,2),mar=c(3,3,2.5,.5),mgp=c(2,1,0))
plot(x,y,ylim=range(y,theta[1,]),col=3)
if(!is.null(u)) lines(x,u,col=2)
lines(x,theta[1,],lwd=2)
title(paste("Reconstruction  h=",signif(hinit,3)))
plot(x,bi[1,],type="l")
title("Sum of weights")
}
if(!is.null(u))
cat("bandwidth: ",signif(hinit,3),"   MSE: ",mean((theta[1,]-u)^2),
      "   MAE: ",mean(abs(theta[1,]-u)),"\n")
if(demo) readline("Press return")
# now run aws-cycle
hakt <- hinit*hincr
while(hakt<=hmax){
z <- .Fortran("lpawsuni",
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
              as.logical(symmetric),PACKAGE="aws")[c("ai","bi","bi0")]
    ai <- (1-eta)*z$ai + eta * ai
    bi <- matrix((1-eta)*z$bi + eta * bi,dp2,n)
    bi0 <- (1-eta)*z$bi0 + eta * bi0
    theta <- matrix(.Fortran("mpawsuni",
                  as.integer(n),
                  as.integer(dp1),
                  as.integer(dp2),
                  as.double(ai),
                  as.double(bi),
                  theta=as.double(theta),
                  double(dp1*dp1),PACKAGE="aws")$theta,dp1,n)
if(graph){
plot(x,y,ylim=range(y,theta[1,]),col=3)
if(!is.null(u)) lines(x,u,col=2)
lines(x,theta[1,],lwd=2)
title(paste("Reconstruction  h=",signif(hakt,3)))
plot(x,bi[1,],type="l")
title("Sum of weights")
}
if(!is.null(u)) 
cat("bandwidth: ",signif(hakt,3),"   MSE: ",mean((theta[1,]-u)^2),
    "   MAE: ",mean(abs(theta[1,]-u)),"\n")
if(demo) readline("Press return")
hakt <- hakt*hincr
gc()
}
}
      if(gridded &&  form=="bi" && p==0){
###                                                                       
###             gridded      bi   p=0                                     
###                                                                       
bi <- ai <- theta <- matrix(0,n1,n2)
if(is.null(hinit)||hinit<1) hinit <- 1
#  first initialize                                                       
z <- .Fortran("iawsbi",
              as.double(y),
              as.integer(n1),
              as.integer(n2),
              as.double(hinit),
              bi=as.double(bi),
              ai=as.double(ai),
              as.double(kernl),
              as.double(wghts),PACKAGE="aws")[c("ai","bi")]
bi <- matrix(z$bi,n1,n2)
ai <- matrix(z$ai,n1,n2)
theta <- ai/bi
if(graph){
par(mfrow=c(1,3),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Observed Image")
image(theta,col=gray((0:255)/255),zlim=range(y),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hinit,3)))
image(bi,col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Sum of weights")
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
z <- .Fortran("lawsbi",
              as.double(y),
              as.integer(n1),
              as.integer(n2),
              as.double(hakt),
              as.double(lamakt),
              as.double(theta),
              bi=as.double(bi),
              ai=as.double(ai),
              as.double(kernl),
              as.double(kerns),
              as.logical(symmetric),
              as.double(wghts),PACKAGE="aws")[c("ai","bi")]
ai <- (1-eta)*z$ai + eta * ai
bi <- matrix((1-eta)*z$bi + eta * bi,n1,n2)
theta  <- matrix(ai / bi, n1, n2)
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Observed Image")
image(theta,col=gray((0:255)/255),zlim=range(y),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)))
image(bi,col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Sum of weights")
if(!is.null(u)) cat("bandwidth: ",signif(hakt,3),"   MSE: ",
                    mean((theta-u)^2),"   MAE: ",mean(abs(theta-u)),"\n")
if(demo) readline("Press return")
hakt <- hakt*hincr
gc()
}
} else
{
#   run all iterations in one call
theta <- .Fortran("gawsbi",
              as.double(y),
              as.integer(n1),
              as.integer(n2),
              as.double(hinit),
              as.double(hincr),
              as.double(hmax),
              as.double(lamakt),
              as.double(eta),
              theta=as.double(theta),
              as.double(bi),
              as.double(ai),
              as.double(kernl),
              as.double(kerns),
              as.double(bi),
              as.logical(symmetric),
              as.double(wghts),PACKAGE="aws")$theta
theta <- matrix(theta,n1,n2)
}
}
      if(gridded &&  form=="bi" && p>0){
###
###             gridded      bi    p=1,2
###
dp1 <- switch(p+1,1,3,6)
dp2 <- switch(p+1,1,6,15)
if(symmetric) dpm <- dp1*(dp1+1)/2 else dpm <- 1
bi <- matrix(0,dp2,n)
theta <- ai <- matrix(0,dp1,n)
ind <- matrix(c(1, 2, 3, 4, 5, 6,
                2, 4, 5, 7, 8, 9,
                3, 5, 6, 8, 9,10,
                4, 7, 8,11,12,13,
                5, 8, 9,12,13,14,
                6, 9,10,13,14,15),6,6)[1:dp1,1:dp1]
if(is.null(hinit)||hinit<p+1.25) hinit <- p+1.25
#  first initialize
z <- .Fortran("ipawsbi",
              as.integer(n1),
              as.integer(n2),
              as.integer(dp1),
              as.integer(dp2),
              as.double(y),
              as.double(hinit),
              bi=as.double(bi),
              ai=as.double(ai),
              theta=as.double(theta),
              as.double(kernl),
              double(dp1*dp1),
              double(dp2),
              double(dp1),
              as.integer(ind),
              as.double(wghts),PACKAGE="aws")[c("ai","bi","theta")]
theta <- array(z$theta,c(dp1,n1,n2))
bi <- bi0 <- array(z$bi,c(dp2,n1,n2))
ai <- z$ai
if(graph){
par(mfrow=c(1,3),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Observed Image")
image(theta[1,,],col=gray((0:255)/255),zlim=range(y),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hinit,3)))
image(bi[1,,],col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Sum of weights")
}
if(!is.null(u))
cat("bandwidth: ",signif(hinit,3),"   MSE: ",mean((theta[1,,]-u)^2),
    "   MAE: ",mean(abs(theta[1,,]-u)),"\n")
if(demo) readline("Press return")
# now run aws-cycle                                                         
hakt <- hinit*hincr
while(hakt<=hmax){
z <- .Fortran("lpawsbi",
              as.integer(n1),
              as.integer(n2),
              as.integer(dp1),
              as.integer(dp2),
              as.double(y),
              as.double(theta),
              as.double(bi),
              bi=as.double(bi),
              bi0=as.double(bi0),
              ai=as.double(ai),
              as.double(lamakt),
              as.double(tau),
              as.double(hakt),
              as.double(kernl),
              as.double(kerns),
              double(dp1*dp1),
              double(dp1*dp1),
              double(dp1*dp1),
              double(dp1),
              double(dp1),
              double(dp2),
              double(dp2),
              double(dp2),
              double(dp1),
              as.logical(symmetric),
              as.integer(ind),
              as.double(wghts),PACKAGE="aws")[c("ai","bi","bi0")]
    ai <- (1-eta)*z$ai + eta * ai
    bi <- array((1-eta)*z$bi + eta * bi,c(dp2,n1,n2))
    bi0 <- (1-eta)*z$bi0 + eta * bi0
    theta <- array(.Fortran("mpawsbi",
                  as.integer(n),
                  as.integer(dp1),
                  as.integer(dp2),
                  as.double(ai),
                  as.double(bi),
                  theta=as.double(theta),
                  double(dp1*dp1),
                  as.integer(ind),PACKAGE="aws")$theta,c(dp1,n1,n2))
if(graph){
par(mfrow=c(1,3),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Observed Image")
image(theta[1,,],col=gray((0:255)/255),zlim=range(y),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)))
image(bi[1,,],col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Sum of weights")
}
if(!is.null(u)) 
   cat("bandwidth: ",signif(hakt,3),"   MSE: ",mean((theta[1,,]-u)^2),
       "   MAE: ",mean(abs(theta[1,,]-u)),"\n")
hakt <- hakt*hincr
gc()
}
}
      if(gridded &&  form=="tri" && p==0){
###                                                                       
###             gridded      tri   p=0                                    
###                                                                       
if(is.null(hinit)||hinit<1) hinit <- 1
#  first initialize
z <- .Fortran("iawstri",
              as.double(y),
              as.integer(n1),
              as.integer(n2),
              as.integer(n3),
              as.double(hinit),
              bi=double(n),
              ai=double(n),
              as.double(kernl),
              as.double(wghts),PACKAGE="aws")[c("ai","bi")]
bi <- array(z$bi,c(n1,n2,n3))
theta <- array(z$ai/bi,c(n1,n2,n3))
if(graph){
par(mfrow=c(1,3),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y[,,1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Observed Image")
image(theta[,,1],col=gray((0:255)/255),zlim=range(y),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hinit,3)))
image(bi[,,1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Sum of weights")
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
z <- .Fortran("lawstri",
              as.double(y),
              as.integer(n1),
              as.integer(n2),
              as.integer(n3),
              as.double(hakt),
              as.double(lamakt),
              as.double(theta),
              bi=as.double(bi),
              ai=double(n),
              as.double(kernl),
              as.double(kerns),
              as.logical(symmetric),
              as.double(wghts),PACKAGE="aws")[c("ai","bi")]
ai <- (1-eta)*z$ai + eta * bi * theta
bi <- array((1-eta)*z$bi + eta * bi,c(n1,n2,n3))
theta  <- array(ai / bi, c(n1,n2,n3))
rm(ai)
gc()
image(y[,,1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Observed Image")
image(theta[,,1],col=gray((0:255)/255),zlim=range(y),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)))
image(bi[,,1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Sum of weights")
if(!is.null(u)) cat("bandwidth: ",signif(hakt,3),"   MSE: ",
                    mean((theta-u)^2),"   MAE: ",mean(abs(theta-u)),"\n")
if(demo) readline("Press return")
hakt <- hakt*hincr
gc()
}
} else
{
#   run all iterations in one call
theta <- array(.Fortran("gawstri",
              as.double(y),
              as.integer(n1),
              as.integer(n2),
              as.integer(n3),
              as.double(hinit),
              as.double(hincr),
              as.double(hmax),
              as.double(lamakt),
              as.double(eta),
              theta=as.double(theta),
              as.double(bi),
              double(n),
              as.double(kernl),
              as.double(kerns),
              as.double(bi),
              as.logical(symmetric),
              as.double(wghts),PACKAGE="aws")$theta, c(n1,n2,n3))
}
}
      if( form=="multi" ){
###
###                        multi (nongridded)    p==0 or p==1
###
dp1 <- 1+p*px
dp2 <- dp1*(dp1+1)/2
bi <- matrix(0,dp2,n)
theta <- ai <- matrix(0,dp1,n)
if(NN){
if(is.null(hinit)||hinit<(p+1)) hinit <- p+1
info <- 1
while(info>0){
ihinit <- trunc(hinit)
z <- .Fortran("ipawsmnn",
              as.integer(n),
              as.integer(px),
              as.integer(dp1),
              as.integer(dp2),
              as.double(x),
              as.double(y),
              as.integer(neighbors[1:ihinit,]),
              as.integer(ihinit),
              as.double(hinit),
              bi=as.double(bi),
              ai=as.double(ai),
              theta=as.double(theta),
              as.double(kernl),
              double(dp1*dp1),
              double(dp1),
              info=as.integer(info),PACKAGE="aws")[c("ai","bi","theta","info")]
info <- z$info
hinit <- hinit+1
}
theta <- matrix(z$theta,dp1,n)
bi <- bi0 <- matrix(z$bi,dp2,n)
ai <- z$ai
if(!is.null(u))
cat("bandwidth: ",signif(hinit,3),"   MSE: ",mean((theta[1,]-u)^2),
                "   MAE: ",mean(abs(theta[1,]-u)),"\n")
# now run aws-cycle
hakt <- hinit*hincr
while(hakt<=hmax){
ihakt <- min(ihmax,trunc(hakt))
z <- .Fortran("lpawsmnn",
              as.integer(n),
              as.integer(px),
              as.integer(dp1),
              as.integer(dp2),
              as.double(x),
              as.double(y),
              as.integer(neighbors[1:ihakt,]),
              as.integer(ihakt),
              as.double(theta),
              as.double(bi),
              bi=as.double(bi),
              bi0=as.double(bi0),
              as.double(ai),
              ai=as.double(ai),
              as.double(lamakt),
              as.double(tau),
              as.double(hakt),
              as.double(kernl),
              as.double(kerns),
              double(dp1*dp1),
              double(dp1*dp1),
              double(dp1*dp1),
              double(dp1),
              double(dp1),
              double(dp1),
              as.logical(symmetric),PACKAGE="aws")[c("ai","bi","bi0")]
    ai <- (1-eta)*z$ai + eta * ai
    bi <- matrix((1-eta)*z$bi + eta * bi,dp2,n)
    bi0 <- (1-eta)*z$bi0 + eta * bi0
    theta <- matrix(.Fortran("mpawsmul",
                             as.integer(n),
                             as.integer(dp1),
                             as.integer(dp2),
                             as.double(ai),
                             as.double(bi),
                             as.double(theta),
                             double(dp1*dp1),PACKAGE="aws")[[6]],dp1,n)
if(!is.null(u))
cat("bandwidth: ",signif(hakt,3),"   MSE: ",mean((theta[1,]-u)^2),
                 "   MAE: ",mean(abs(theta[1,]-u)),"\n")
hakt <- hakt*hincr
gc()
}
} else {
dpd <- dp1+1
if(is.null(hinit)) hinit <- maxdist[dpd]
info <- 1
while(info>0){
if(hinit<=maxdist[dpd]) hinit <- maxdist[dpd]
ihinit <- sum(maxdist<=hinit)
z <- .Fortran("ipawsmul",
              as.integer(n),
              as.integer(px),
              as.integer(dp1),
              as.integer(dp2),
              as.double(x),
              as.double(y),
              as.integer(neighbors[1:ihinit,]),
              as.double(distmat[1:ihinit,]),
              as.integer(ihinit),
              as.double(hinit),
              bi=as.double(bi),
              ai=as.double(ai),
              theta=as.double(theta),
              as.double(kernl),
              double(dp1*dp1),
              double(dp1),
              info=integer(1),PACKAGE="aws")[c("ai","bi","theta","info")]
info <- z$info
dpd <- dpd+1
}
theta <- matrix(z$theta,dp1,n)
bi <- bi0 <- matrix(z$bi,dp2,n)
ai <- z$ai
if(!is.null(u))
cat("bandwidth: ",signif(hinit,3),"   MSE: ",mean((theta[1,]-u)^2),
                "   MAE: ",mean(abs(theta[1,]-u)),"\n")
# now run aws-cycle
hakt <- hinit*hincr
while(hakt<=hmax){
ihakt <- sum(maxdist<=hakt)
z <- .Fortran("lpawsmul",
              as.integer(n),
              as.integer(px),
              as.integer(dp1),
              as.integer(dp2),
              as.double(x),
              as.double(y),
              as.integer(neighbors[1:ihakt,]),
              as.double(distmat[1:ihakt,]),
              as.integer(ihakt),
              as.double(theta),
              as.double(bi),
              bi=as.double(bi),
              bi0=as.double(bi0),
              as.double(ai),
              ai=as.double(ai),
              as.double(lamakt),
              as.double(tau),
              as.double(hakt),
              as.double(kernl),
              as.double(kerns),
              double(dp1*dp1),
              double(dp1*dp1),
              double(dp1*dp1),
              double(dp1),
              double(dp1),
              double(dp1),
              as.logical(symmetric),PACKAGE="aws")[c("ai","bi","bi0")]
    ai <- (1-eta)*z$ai + eta * ai
    bi <- matrix((1-eta)*z$bi + eta * bi,dp2,n)
    bi0 <- (1-eta)*z$bi0 + eta * bi0
    theta <- matrix(.Fortran("mpawsmul",
                             as.integer(n),
                             as.integer(dp1),
                             as.integer(dp2),
                             as.double(ai),
                             as.double(bi),
                             as.double(theta),
                             double(dp1*dp1),PACKAGE="aws")[[6]],dp1,n)
if(!is.null(u))
cat("bandwidth: ",signif(hakt,3),"   MSE: ",mean((theta[1,]-u)^2),
                  "   MAE: ",mean(abs(theta[1,]-u)),"\n")
hakt <- hakt*hincr
gc()
}
}
}
###
###            end cases
###
z<-list(theta=theta,y=y,x=x,call=args)
class(z)<-"aws"
z
}
