#
#    R - function  laws  for likelihood  based  Adaptive Weights Smoothing (AWS)
#    for local constant Gaussian, Bernoulli, Exponential, Poisson, Weibull and  
#    Volatility models                                                          
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
laws <- function(y,x=NULL,qlambda=NULL,eta=0.5,lkern="Triangle",model="Poisson",
                 shape=NULL,hinit=NULL,hincr=NULL,hmax=10,NN=FALSE,u=NULL,
                 graph=FALSE,demo=FALSE,symmetric=TRUE)
{ 
#
#    first check arguments and initialize                                 
#
args <- match.call()
eps <- 1.e-10
if(is.null(qlambda)) if(symmetric) qlambda <- switch(model,
                                                     Gaussian=.985,
                                                     Bernoulli=.972,
                                                     Exponential=.972,
                                                     Poisson=.980,
                                                     Weibull=.972,
                                                     Volatility=.972)
                     else qlambda <- switch(model,
                                                     Gaussian=.966,
                                                     Bernoulli=.953,
                                                     Exponential=.914,
                                                     Poisson=.958,
                                                     Weibull=.914,
                                                     Volatility=.914)
if(qlambda>=1 || qlambda<.6) return("Inappropriate value of qlambda")
if(eta<eps || eta>=1) return("Inappropriate value of eta")
if(model!="Gaussian"&&model!="Bernoulli"&&model!="Exponential"&&
   model!="Poisson"&&model!="Weibull"&&model!="Volatility") 
   return(paste("specified model ",model," not yet implemented"))
#
#    generate kernel on a grid  and set lambda
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
#      get lambda as quantile of appropriate chisq,                       
#                rescale to be consistent with the paper in  lamakt       
#
lamakt <- 5*qchisq(qlambda,1)
if(model=="Gaussian") lamakt <- lamakt*shape*2
#
#   specify which statistics are needed and transform data if necessary   
#
logtheta <- switch(model,Gaussian=FALSE,Bernoulli=TRUE,Exponential=TRUE,
                   Poisson=TRUE,Weibull=TRUE,Volatility=TRUE)
logctheta <- switch(model,Gaussian=FALSE,Bernoulli=TRUE,Exponential=FALSE,
                    Poisson=FALSE,Weibull=FALSE,Volatility=FALSE)
if(model=="Weibull" && (is.null(shape) || shape<=0)) 
   return("Shape parameter for Weibull has to be positive")
if(model=="Gaussian" && (is.null(shape) || shape<=0)) 
   return("Variance (shape) for Gaussian errors has to be positive")
weibull <- FALSE
shape <- 1
if(model=="Weibull") {
model <- "Exponential"
y <- y^shape
weibull <- TRUE
}
if(model=="Volatility"){
model <- "Exponential"
y <- y^2
lamakt <- 2*lamakt 
# this accounts for the additional 1/2 in Q(\hat{theta},theta)
weibull <- TRUE
shape <- 2
}
if(demo&& !graph) graph <- TRUE
# now check which procedure is appropriate
gridded <- is.null(x) 
if(gridded){
##  this is the version on a grid
dy <- dim(y)
if(is.null(dy)) {
   form <- "uni"
   ddim  <- 1
   n <- length(y)
}
if(length(dy)==2){
   form <- "bi"
   ddim  <- 2
n1 <- dy[1]
n2 <- dy[2]
n <- n1*n2
}
if(length(dy)==3){
   form <- "tri"
   ddim  <- 3
n1 <- dy[1]
n2 <- dy[2]
n3 <- dy[3]
n <- n1*n2*n3
}
if(length(dy)>3) 
   return("AWS for more than 3 dimensional grids is not implemented")
} else { 
# not gridded
dx <- dim(x)   
ddim <- 1
if(is.null(dx)&&NN) {
#
#    order data by order of x
#
    form <- "uni"
    n <- length(x) 
    if(n!=length(y)) return("incompatible lengths of x and y")
    ox <- order(x)
    x <- x[ox]
    y <- y[ox]
}else {
   if(is.null(dx)){
      px <- 1
      n <- length(x)
   }else{
   px <- dx[1]
   n <- dx[2]
   }
   form <- "multi"
   if(n!=length(y)) return("incompatible dimensions of x and y")
   weights <- rep(1,px)
#
#  now generate matrix of nearest neighbors
#  hmax is interpreted as maximal number of neighbors
#
   if(NN){
   ihmax <- trunc(hmax)
   if(ihmax>n) ihmax <- n
   neighbors <- matrix(0,ihmax,n)
   for (i in 1:n) {
      adist <- weights%*%((x-x[,i])^2)
      neighbors[,i] <- order(adist)[1:ihmax]
      }
   } else {
   ihmax <- n
   ddim <- px
   neighbors <- distmat <- matrix(0,n,n)
   for (i in 1:n) {
      adist <- weights%*%((x-x[,i])^2)
      od <- order(adist)
      distmat[,i] <- adist[od]
      neighbors[,i] <- od
      }
#  now reduce memory used to whats needed
   gc()
   distmat <- sqrt(distmat)
   maxdist <- apply(distmat,1,max)
   meandist <- apply(distmat,1,mean)
   mindist <- apply(distmat,1,min)
   ihmax <- sum(mindist<=hmax)
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
#     now set hincr if not provided                                       
#
if(is.null(hincr)) hincr <- 1.25^(1/ddim)  
#
#    get a global estimate if they are needed for regularization          
#
if(logtheta) gtheta <- mean(y)
#
#    now select the correct aws-procedure                                 
#
#   cases:    gridded      uni                                            
#             gridded      bi                                             
#             gridded      tri                                            
#             !gridded     multi                                          
#             !gridded     multi, Nearest Neighbor                        
#
if(gridded &&  form=="uni" ){
###                                                                       
###              gridded     uni                                          
###                                                                       
###     this should run a little faster than the nongridded version       
###                                                                       
bi <- ai <- theta <- numeric(n)
if(is.null(hinit)||hinit<1) hinit <- 1
#  first initialize
z <- .Fortran("iawsuni",
              as.double(y),
              as.integer(n),
              as.double(hinit),
              bi=as.double(bi),
              ai=as.double(ai),
              as.double(kernl))[c("bi","ai")]
bi <- z$bi
ai <- z$ai
if(logtheta) {
bi <- (1-eta)*bi+eta
ai <- (1-eta)*ai+eta*gtheta
}
theta <- ai/bi
if(logtheta) ltheta <- log(theta+eps*max(theta))
if(logctheta) lctheta <- log(1.e0-theta+eps*max(1.e0-theta))
if(weibull) theta <- theta^(1/shape)
if(graph){
par(mfrow=c(1,2),mar=c(3,3,2.5,.5),mgp=c(2,1,0))
plot(y^(1/shape),ylim=range(y^(1/shape),theta),col=3)
if(!is.null(u)) lines(u,col=2)
lines(theta,lwd=2)
title(paste("Reconstruction  h=",signif(hinit,3)))
plot(bi,type="l")
title("Sum of weights")
}
if(!is.null(u)) cat("bandwidth: ",signif(hinit,3),"   MSE: ",
                    mean((theta-u)^2),"   MAE: ",mean(abs(theta-u)),"\n")
if(demo) readline("Press return")
if(weibull) theta <- theta^(shape)
# now run aws-cycle                                                       
hakt <- hinit*hincr
if(graph){
#
#   run single steps to display intermediate results                      
#
while(hakt<=hmax){
z <- switch(model,
            Gaussian=.Fortran("lawsuni",
                              as.double(y),
                              as.integer(n),
                              as.double(hakt),
                              as.double(lamakt),
                              as.double(theta),
                              bi=as.double(bi),
                              ai=as.double(ai),
                              as.double(kernl),
                              as.double(kerns),
                              as.logical(symmetric))[c("bi","ai")],
            Bernoulli=.Fortran("lberuni",
                               as.double(y),
                               as.integer(n),
                               as.double(hakt),
                               as.double(lamakt),
                               as.double(theta),
                               as.double(ltheta),
                               as.double(lctheta),
                               bi=as.double(bi),
                               ai=as.double(ai),
                               as.double(kernl),
                               as.double(kerns),
                              as.logical(symmetric))[c("bi","ai")],
            Poisson=.Fortran("lpoiuni",
                             as.double(y),
                             as.integer(n),
                             as.double(hakt),
                             as.double(lamakt),
                             as.double(theta),
                             as.double(ltheta),
                             bi=as.double(bi),
                             ai=as.double(ai),
                             as.double(kernl),
                             as.double(kerns),
                              as.logical(symmetric))[c("bi","ai")],
            Exponential=.Fortran("lexpuni",
                             as.double(y),
                             as.integer(n),
                             as.double(hakt),
                             as.double(lamakt),
                             as.double(theta),
                             as.double(ltheta),
                             bi=as.double(bi),
                             ai=as.double(ai),
                             as.double(kernl),
                             as.double(kerns),
                              as.logical(symmetric))[c("bi","ai")])
ai <- (1-eta)*z$ai + eta * ai
bi <- (1-eta)*z$bi + eta * bi
theta  <- ai / bi
if(logtheta) ltheta <- log(theta+eps*max(theta))
if(logctheta) lctheta <- log(1.e0-theta+eps*max(1.e0-theta))
if(weibull) theta <- theta^(1/shape)
plot(y^(1/shape),ylim=range(y^(1/shape),theta),col=3)
if(!is.null(u)) lines(u,col=2)
lines(theta,lwd=2)
title(paste("Reconstruction  h=",signif(hakt,3)))
plot(bi,type="l")
title("Sum of weights")
if(!is.null(u)) cat("bandwidth: ",signif(hakt,3),"   MSE: ",
                    mean((theta-u)^2),"   MAE: ",mean(abs(theta-u)),"\n")
if(demo) readline("Press return")
if(weibull) theta <- theta^(shape)
hakt <- hakt*hincr
gc()
}
if(weibull) theta <- theta^(1/shape)
} else
{
#   run all iterations in one call                                        
theta <- switch(model,
            Gaussian=.Fortran("gawsuni",
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
                              as.double(ai),
                              as.logical(symmetric))$theta,
            Bernoulli=.Fortran("gberuni",
                              as.double(y),
                              as.integer(n),
                              as.double(hinit),
                              as.double(hincr),
                              as.double(hmax),
                              as.double(lamakt),
                              as.double(eta),
                              theta=as.double(theta),
                              as.double(ltheta),
                              as.double(lctheta),
                              as.double(bi),
                              as.double(ai),
                              as.double(kernl),
                              as.double(kerns),
                              as.double(bi),
                              as.double(ai),
                              as.logical(symmetric))$theta,
            Poisson= .Fortran("gpoiuni",
                              as.double(y),
                              as.integer(n),
                              as.double(hinit),
                              as.double(hincr),
                              as.double(hmax),
                              as.double(lamakt),
                              as.double(eta),
                              theta=as.double(theta),
                              as.double(ltheta),
                              as.double(bi),
                              as.double(ai),
                              as.double(kernl),
                              as.double(kerns),
                              as.double(bi),
                              as.double(ai),
                              as.logical(symmetric))$theta,
            Exponential=.Fortran("gexpuni",
                              as.double(y),
                              as.integer(n),
                              as.double(hinit),
                              as.double(hincr),
                              as.double(hmax),
                              as.double(lamakt),
                              as.double(eta),
                              theta=as.double(theta),
                              as.double(ltheta),
                              as.double(bi),
                              as.double(ai),
                              as.double(kernl),
                              as.double(kerns),
                              as.double(bi),
                              as.double(ai),
                              as.logical(symmetric))$theta)
if(weibull) theta <- theta^(1/shape)
}
}
      if(gridded &&  form=="bi" ){
###                                                                       
###             gridded      bi                                           
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
              as.double(kernl))[c("bi","ai")]
bi <- z$bi
ai <- z$ai
if(logtheta) {
bi <- (1-eta)*bi+eta
ai <- (1-eta)*ai+eta*gtheta
}
theta <- matrix(ai/bi,n1,n2)
if(logtheta) ltheta <- log(theta+eps*max(theta))
if(logctheta) lctheta <- log(1.e0-theta+eps*max(1.e0-theta))
bi <- matrix(bi,n1,n2)
if(weibull) theta <- theta^(1/shape)
if(graph){
par(mfrow=c(1,3),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Observed Image")
image(theta,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hinit,3)))
image(bi,col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Sum of weights")
}
if(!is.null(u)) cat("bandwidth: ",signif(hinit,3),"   MSE: ",
                    mean((theta-u)^2),"   MAE: ",mean(abs(theta-u)),"\n")
if(demo) readline("Press return")
if(weibull) theta <- theta^(shape)
# now run aws-cycle                                                       
hakt <- hinit*hincr
if(graph){
#
#   run single steps to display intermediate results                      
#
while(hakt<=hmax){
z <- switch(model,
            Gaussian=.Fortran("lawsbi",
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
                              as.logical(symmetric))[c("bi","ai")],
            Bernoulli=.Fortran("lberbi",
                              as.double(y),
                              as.integer(n1),
                              as.integer(n2),
                              as.double(hakt),
                              as.double(lamakt),
                              as.double(theta),
                              as.double(ltheta),
                              as.double(lctheta),
                              bi=as.double(bi),
                              ai=as.double(ai),
                              as.double(kernl),
                              as.double(kerns),
                              as.logical(symmetric))[c("bi","ai")],
            Poisson=.Fortran("lpoibi",
                              as.double(y),
                              as.integer(n1),
                              as.integer(n2),
                              as.double(hakt),
                              as.double(lamakt),
                              as.double(theta),
                              as.double(ltheta),
                              bi=as.double(bi),
                              ai=as.double(ai),
                              as.double(kernl),
                              as.double(kerns),
                              as.logical(symmetric))[c("bi","ai")],
            Exponential=.Fortran("lexpbi",
                              as.double(y),
                              as.integer(n1),
                              as.integer(n2),
                              as.double(hakt),
                              as.double(lamakt),
                              as.double(theta),
                              as.double(ltheta),
                              bi=as.double(bi),
                              ai=as.double(ai),
                              as.double(kernl),
                              as.double(kerns),
                              as.logical(symmetric))[c("bi","ai")])
ai <- (1-eta)*z$ai + eta * ai
bi <- matrix((1-eta)*z$bi + eta * bi,n1,n2)
theta  <- matrix(ai / bi, n1, n2)
if(logtheta) ltheta <- log(theta+eps*max(theta))
if(logctheta) lctheta <- log(1.e0-theta+eps*max(1.e0-theta))
if(weibull) theta <- theta^(1/shape)
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Observed Image")
image(theta,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)))
image(bi,col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Sum of weights")
if(!is.null(u)) cat("bandwidth: ",signif(hakt,3),"   MSE: ",
                    mean((theta-u)^2),"   MAE: ",mean(abs(theta-u)),"\n")
if(demo) readline("Press return")
hakt <- hakt*hincr
if(weibull) theta <- theta^(shape)
gc()
}
if(weibull) theta <- theta^(1/shape)
} else
{
#   run all iterations in one call                                        
theta <- switch(model,
            Gaussian=.Fortran("gawsbi",
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
                              as.double(ai),
                              as.logical(symmetric))$theta,
            Bernoulli=.Fortran("gberbi",
                              as.double(y),
                              as.integer(n1),
                              as.integer(n2),
                              as.double(hinit),
                              as.double(hincr),
                              as.double(hmax),
                              as.double(lamakt),
                              as.double(eta),
                              theta=as.double(theta),
                              as.double(ltheta),
                              as.double(lctheta),
                              as.double(bi),
                              as.double(ai),
                              as.double(kernl),
                              as.double(kerns),
                              as.double(bi),
                              as.double(ai),
                              as.logical(symmetric))$theta,
            Poisson=.Fortran("gpoibi",
                              as.double(y),
                              as.integer(n1),
                              as.integer(n2),
                              as.double(hinit),
                              as.double(hincr),
                              as.double(hmax),
                              as.double(lamakt),
                              as.double(eta),
                              theta=as.double(theta),
                              as.double(ltheta),
                              as.double(bi),
                              as.double(ai),
                              as.double(kernl),
                              as.double(kerns),
                              as.double(bi),
                              as.double(ai),
                              as.logical(symmetric))$theta,
            Exponential=.Fortran("gexpbi",
                              as.double(y),
                              as.integer(n1),
                              as.integer(n2),
                              as.double(hinit),
                              as.double(hincr),
                              as.double(hmax),
                              as.double(lamakt),
                              as.double(eta),
                              theta=as.double(theta),
                              as.double(ltheta),
                              as.double(bi),
                              as.double(ai),
                              as.double(kernl),
                              as.double(kerns),
                              as.double(bi),
                              as.double(ai),
                              as.logical(symmetric))$theta)
theta <- matrix(theta,n1,n2)
if(weibull) theta <- theta^(1/shape)
}
}
      if(gridded &&  form=="tri" ){
###                                                                       
###             gridded      tri                                          
###                                                                       
bi <- ai <- theta <- array(0,c(n1,n2,n3))
if(is.null(hinit)||hinit<1) hinit <- 1
#  first initialize
z <- .Fortran("iawstri",
              as.double(y),
              as.integer(n1),
              as.integer(n2),
              as.integer(n3),
              as.double(hinit),
              bi=as.double(bi),
              ai=as.double(ai),
              as.double(kernl))[c("bi","ai")]
bi <- z$bi
ai <- z$ai
if(logtheta) {
bi <- (1-eta)*bi+eta
ai <- (1-eta)*ai+eta*gtheta
}
theta <- array(ai/bi,c(n1,n2,n3))
if(logtheta) ltheta <- log(theta+eps*max(theta))
if(logctheta) lctheta <- log(1.e0-theta+eps*max(1.e0-theta))
bi <- array(z$bi,c(n1,n2,n3))
if(weibull) theta <- theta^(1/shape)
if(graph){
par(mfrow=c(1,3),mar=c(1,1,3,.25),mgp=c(2,1,0))
image(y[,,1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Observed Image")
image(theta[,,1],col=gray((0:255)/255),zlim=range(y),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hinit,3)))
image(bi[,,1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Sum of weights")
}
if(!is.null(u)) cat("bandwidth: ",signif(hinit,3),"   MSE: ",
                    mean((theta-u)^2),"   MAE: ",mean(abs(theta-u)),"\n")
if(demo) readline("Press return")
if(weibull) theta <- theta^(shape)
# now run aws-cycle                                                       
hakt <- hinit*hincr
if(graph){
#
#   run single steps to display intermediate results                      
#
while(hakt<=hmax){
z <- switch(model,
            Gaussian=.Fortran("lawstri",
                               as.double(y),
                               as.integer(n1),
                               as.integer(n2),
                               as.integer(n3),
                               as.double(hakt),
                               as.double(lamakt),
                               as.double(theta),
                               bi=as.double(bi),
                               ai=as.double(ai),
                               as.double(kernl),
                               as.double(kerns),
                              as.logical(symmetric))[c("bi","ai")],
            Bernoulli=.Fortran("lbertri",
                               as.double(y),
                               as.integer(n1),
                               as.integer(n2),
                               as.integer(n3),
                               as.double(hakt),
                               as.double(lamakt),
                               as.double(theta),
                               as.double(ltheta),
                               as.double(lctheta),
                               bi=as.double(bi),
                               ai=as.double(ai),
                               as.double(kernl),
                               as.double(kerns),
                              as.logical(symmetric))[c("bi","ai")],
            Poisson=.Fortran("lpoitri",
                              as.double(y),
                              as.integer(n1),
                              as.integer(n2),
                              as.integer(n3),
                              as.double(hakt),
                              as.double(lamakt),
                              as.double(theta),
                              as.double(ltheta),
                              bi=as.double(bi),
                              ai=as.double(ai),
                              as.double(kernl),
                              as.double(kerns),
                              as.logical(symmetric))[c("bi","ai")],
            Exponential=.Fortran("lexptri",
                              as.double(y),
                              as.integer(n1),
                              as.integer(n2),
                              as.integer(n3),
                              as.double(hakt),
                              as.double(lamakt),
                              as.double(theta),
                              as.double(ltheta),
                              bi=as.double(bi),
                              ai=as.double(ai),
                              as.double(kernl),
                              as.double(kerns),
                              as.logical(symmetric))[c("bi","ai")])
ai <- (1-eta)*z$ai + eta * bi * theta
bi <- array((1-eta)*z$bi + eta * bi,c(n1,n2,n3))
theta  <- array(ai / bi, c(n1,n2,n3))
if(logtheta) ltheta <- log(theta+eps*max(theta))
if(logctheta) lctheta <- log(1.e0-theta+eps*max(1.e0-theta))
if(weibull) theta <- theta^(1/shape)
image(y[,,1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Observed Image")
image(theta[,,1],col=gray((0:255)/255),zlim=range(y),xaxt="n",yaxt="n")
title(paste("Reconstruction  h=",signif(hakt,3)))
image(bi[,,1],col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Sum of weights")
if(!is.null(u)) cat("bandwidth: ",signif(hakt,3),"   MSE: ",
                    mean((theta-u)^2),"   MAE: ",mean(abs(theta-u)),"\n")
if(demo) readline("Press return")
if(weibull) theta <- theta^(shape)
hakt <- hakt*hincr
gc()
}
if(weibull) theta <- theta^(1/shape)
} else
{
#   run all iterations in one call                                        
theta <- switch(model,
            Gaussian=.Fortran("gawstri",
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
                               as.double(ai),
                               as.double(kernl),
                               as.double(kerns),
                               as.double(bi),
                               as.double(ai),
                              as.logical(symmetric))$theta,
            Bernoulli=.Fortran("gbertri",
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
                               as.double(ltheta),
                               as.double(lctheta),
                               as.double(bi),
                               as.double(ai),
                               as.double(kernl),
                               as.double(kerns),
                               as.double(bi),
                               as.double(ai),
                              as.logical(symmetric))$theta,
            Poisson=.Fortran("gpoitri",
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
                              as.double(ltheta),
                              as.double(bi),
                              as.double(ai),
                              as.double(kernl),
                              as.double(kerns),
                              as.double(bi),
                              as.double(ai),
                              as.logical(symmetric))$theta,
            Exponential=.Fortran("gexptri",
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
                              as.double(ltheta),
                              as.double(bi),
                              as.double(ai),
                              as.double(kernl),
                              as.double(kerns),
                              as.double(bi),
                              as.double(ai),
                              as.logical(symmetric))$theta)
theta <- array(theta, c(n1,n2,n3))
if(weibull) theta <- theta^(1/shape)
}
}
      if( form=="multi" ){
###                                                                       
###                        multi (nongridded)    p==0 or p==1             
###                                                                       
bi <- numeric(n)
theta <- ai <- numeric(n)
if(NN){
if(is.null(hinit)||hinit<2) hinit <- 2
ihinit <- trunc(hinit)
z <- .Fortran("iawsmnn",
              as.integer(n),
              as.double(y),
              as.integer(neighbors[1:ihinit,]),
              as.integer(ihinit),
              as.double(hinit),
              bi=as.double(bi),
              ai=as.double(ai),
              as.double(kernl))[c("bi","ai")]
bi <- z$bi
ai <- z$ai
if(logtheta) {
bi <- (1-eta)*bi+eta
ai <- (1-eta)*ai+eta*gtheta
}
theta <- ai/bi
if(logtheta) ltheta <- log(theta+eps*max(theta))
if(logctheta) lctheta <- log(1.e0-theta+eps*max(1.e0-theta))
if(!is.null(u)){ 
if(weibull) theta <- theta^(1/shape)
cat("bandwidth: ",signif(hinit,3),"   MSE: ",
    mean((theta-u)^2),"   MAE: ",mean(abs(theta-u)),"\n")
if(weibull) theta <- theta^(shape)
}
# now run aws-cycle
hakt <- hinit*hincr
while(hakt<=hmax){
ihakt <- min(ihmax,trunc(hakt))
z <- switch(model,
            Gaussian=.Fortran("lawsmnn",
                               as.integer(n),
                               as.double(y),
                               as.integer(neighbors[1:ihakt,]),
                               as.integer(ihakt),
                               as.double(theta),
                               as.double(bi),
                               bi=as.double(bi),
                               as.double(ai),
                               ai=as.double(ai),
                               as.double(lamakt),
                               as.double(hakt),
                               as.double(kernl),
                               as.double(kerns),
                              as.logical(symmetric))[c("bi","ai")],
              Bernoulli=.Fortran("lbermnn",
                               as.integer(n),
                               as.double(y),
                               as.integer(neighbors[1:ihakt,]),
                               as.integer(ihakt),
                               as.double(theta),
                               as.double(ltheta),
                               as.double(lctheta),
                               as.double(bi),
                               bi=as.double(bi),
                               as.double(ai),
                               ai=as.double(ai),
                               as.double(lamakt),
                               as.double(hakt),
                               as.double(kernl),
                               as.double(kerns),
                              as.logical(symmetric))[c("bi","ai")],
              Poisson=.Fortran("lpoimnn",
                               as.integer(n),
                               as.double(y),
                               as.integer(neighbors[1:ihakt,]),
                               as.integer(ihakt),
                               as.double(theta),
                               as.double(ltheta),
                               as.double(bi),
                               bi=as.double(bi),
                               as.double(ai),
                               ai=as.double(ai),
                               as.double(lamakt),
                               as.double(hakt),
                               as.double(kernl),
                               as.double(kerns),
                              as.logical(symmetric))[c("bi","ai")],
              Exponential=.Fortran("lexpmnn",
                               as.integer(n),
                               as.double(y),
                               as.integer(neighbors[1:ihakt,]),
                               as.integer(ihakt),
                               as.double(theta),
                               as.double(ltheta),
                               as.double(bi),
                               bi=as.double(bi),
                               as.double(ai),
                               ai=as.double(ai),
                               as.double(lamakt),
                               as.double(hakt),
                               as.double(kernl),
                               as.double(kerns),
                              as.logical(symmetric))[c("bi","ai")])
    ai <- (1-eta)*z$ai + eta * ai
    bi <- (1-eta)*z$bi + eta * bi
    theta <- ai/bi
if(logtheta) ltheta <- log(theta+eps*max(theta))
if(logctheta) lctheta <- log(1.e0-theta+eps*max(1.e0-theta))
if(!is.null(u)) {
if(weibull) theta <- theta^(1/shape)
cat("bandwidth: ",signif(hakt,3),"   MSE: ",
    mean((theta-u)^2),"   MAE: ",mean(abs(theta-u)),"\n")
if(weibull) theta <- theta^(shape)
}
hakt <- hakt*hincr
gc()
}
if(weibull) theta <- theta^(1/shape)
} else {
dpd <- 2
if(is.null(hinit)) hinit <- maxdist[dpd]
if(hinit<=meandist[dpd]) hinit <- meandist[dpd]
ihinit <- sum(mindist<=hinit)
z <- .Fortran("iawsmul",
              as.integer(n),
              as.double(y),
              as.integer(neighbors[1:ihinit,]),
              as.double(distmat[1:ihinit,]),
              as.integer(ihinit),
              as.double(hinit),
              bi=as.double(bi),
              ai=as.double(ai),
              as.double(kernl))[c("bi","ai")]
bi <- z$bi
ai <- z$ai
if(logtheta) {
bi <- (1-eta)*bi+eta
ai <- (1-eta)*ai+eta*gtheta
}
theta <- ai/bi
if(logtheta) ltheta <- log(theta+eps*max(theta))
if(logctheta) lctheta <- log(1.e0-theta+eps*max(1.e0-theta))
if(!is.null(u)) {
if(weibull) theta <- theta^(1/shape)
cat("bandwidth: ",signif(hinit,3),"   MSE: ",mean((theta-u)^2),
    "   MAE: ",mean(abs(theta-u)),"\n")
if(weibull) theta <- theta^(shape)
}
# now run aws-cycle
hakt <- hinit*hincr
while(hakt<=hmax){
ihakt <- sum(maxdist<=hakt)
z <- switch(model,
            Gaussian=.Fortran("lawsmul",
                               as.integer(n),
                               as.double(y),
                               as.integer(neighbors[1:ihakt,]),
                               as.double(distmat[1:ihakt,]),
                               as.integer(ihakt),
                               as.double(theta),
                               as.double(bi),
                               bi=as.double(bi),
                               as.double(ai),
                               ai=as.double(ai),
                               as.double(lamakt),
                               as.double(hakt),
                               as.double(kernl),
                               as.double(kerns),
                              as.logical(symmetric))[c("bi","ai")],
              Bernoulli=.Fortran("lbermul",
                               as.integer(n),
                               as.double(y),
                               as.integer(neighbors[1:ihakt,]),
                               as.double(distmat[1:ihakt,]),
                               as.integer(ihakt),
                               as.double(theta),
                               as.double(ltheta),
                               as.double(lctheta),
                               as.double(bi),
                               bi=as.double(bi),
                               as.double(ai),
                               ai=as.double(ai),
                               as.double(lamakt),
                               as.double(hakt),
                               as.double(kernl),
                               as.double(kerns),
                              as.logical(symmetric))[c("bi","ai")],
              Poisson=.Fortran("lpoimul",
                               as.integer(n),
                               as.double(y),
                               as.integer(neighbors[1:ihakt,]),
                               as.double(distmat[1:ihakt,]),
                               as.integer(ihakt),
                               as.double(theta),
                               as.double(ltheta),
                               as.double(bi),
                               bi=as.double(bi),
                               as.double(ai),
                               ai=as.double(ai),
                               as.double(lamakt),
                               as.double(hakt),
                               as.double(kernl),
                               as.double(kerns),
                              as.logical(symmetric))[c("bi","ai")],
              Exponential=.Fortran("lexpmul",
                               as.integer(n),
                               as.double(y),
                               as.integer(neighbors[1:ihakt,]),
                               as.double(distmat[1:ihakt,]),
                               as.integer(ihakt),
                               as.double(theta),
                               as.double(ltheta),
                               as.double(bi),
                               bi=as.double(bi),
                               as.double(ai),
                               ai=as.double(ai),
                               as.double(lamakt),
                               as.double(hakt),
                               as.double(kernl),
                               as.double(kerns),
                              as.logical(symmetric))[c("bi","ai")])
    ai <- (1-eta)*z$ai + eta * ai
    bi <- (1-eta)*z$bi + eta * bi
    theta <- ai/bi
if(logtheta) ltheta <- log(theta+eps*max(theta))
if(logctheta) lctheta <- log(1.e0-theta+eps*max(1.e0-theta))
if(!is.null(u)) {
if(weibull) theta <- theta^(1/shape)
cat("bandwidth: ",signif(hakt,3),"   MSE: ",
    mean((theta-u)^2),"   MAE: ",mean(abs(theta-u)),"\n")
if(weibull) theta <- theta^(shape)
}
hakt <- hakt*hincr
gc()
}
if(weibull) theta <- theta^(1/shape)
}
}
###                                                                       
###            end cases                                                  
###                                                                       
list(theta=theta,y=y,x=x,call=args)
}
