#
#    R - function  awstindex  for tail-index estimation                                                          
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
awstindex <- function(y,qlambda=NULL,eta=0.5,lkern="Triangle",hinit=1,hincr=1.25,hmax=1000,
                 graph=FALSE,symmetric=FALSE){
args <- match.call()
if(!is.null(dim(y))) 
   return("only univariate tail index estimation is implemented")
n <- length(y)
y <- sort(y)[n:1]
x <- (1:(n-1))*log(y[-n]/y[-1])
theta <- laws(x,qlambda=qlambda,model="Exponential",hinit=hinit,hincr=hincr,
     hmax=hmax,graph=graph,symmetric=symmetric)$theta
list(tindex=theta[1],intensity=theta,y=y,call=args)
}
#
#    R - function  awsdens  for local constant density estimation  in 1D, 2D and 3D                                                        
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
awsdens <- function(y,ngrid=NULL,nempty=NULL,qlambda=NULL,eta=0.5,lkern="Triangle",fu=NULL,
                    hinit=1,hincr=1.2,hmax=NULL,graph=FALSE,demo=FALSE,symmetric=TRUE){
# first generate the grid of bin's
args <- match.call()
dy <- dim(y)
u <- NULL
if(is.null(dy)){
   n <- length(y)
   if(is.null(ngrid)) ngrid <- as.integer(2*n)
   if(is.null(nempty)) nempty <- trunc(.1*ngrid)
   dd <- 1
   ry <- range(y)
   dry <- diff(ry)
   ry[1] <- ry[1]-dry/ngrid*nempty
   ry[2] <- ry[2]+dry/ngrid*nempty
   dry <- dry*(1+2/ngrid*nempty)
   bin <- numeric(ngrid)
   ind <- trunc((y-ry[1])/dry*ngrid)+1
   bin[as.integer(levels(factor(ind)))] <- table(ind)
   xgrid <- list(seq(ry[1]+dry/2/ngrid,ry[2]+dry/2/ngrid,length=ngrid))
   if(!is.null(fu)) u <- fu(xgrid)
   } else {
# assume the rows to contain components
   dd <- dy[1] 
   n <- dy[2]
   if(is.null(ngrid)) ngrid <- as.integer(2*n^(1/dd))
   if(is.null(nempty)) nempty <- trunc(.1*ngrid)
   ry <- apply(y,1,range)
   dry <- ry[2,]-ry[1,]
   ry[1,] <- ry[1,]-dry/ngrid*nempty
   ry[2,] <- ry[2,]+dry/ngrid*nempty
   dry <- dry*(1+2/ngrid*nempty)
   if(length(ngrid)==1) ngrid <- rep(ngrid,dd)
   if(length(ngrid)!=dd) return("incompatible length of ngrid")
   if(dd>3) return("not implemented for more than three dimensions")
   if(dd==2) bin <- matrix(0,ngrid[1],ngrid[2])
   if(dd>2) bin <- array(0,ngrid)
   ind <- matrix(0,n,dd)
   for(i in 1:dd) 
       ind[,i] <- trunc((y[i,]-ry[1,i])/dry[i]*ngrid[i])+1
   if(dd==2) {
       bin[as.integer(levels(factor(ind[,1]))),
       as.integer(levels(factor(ind[,2])))] <- table(ind[,1],ind[,2]) 
       xgrid <- list(seq(ry[1,1]+dry[1]/2/ngrid[1],ry[2,1]+dry[1]/2/ngrid[1],length=ngrid[1]),
                      seq(ry[1,2]+dry[2]/2/ngrid[2],ry[2,2]+dry[2]/2/ngrid[2],length=ngrid[2]))
       if(!is.null(fu)) u <- fu(xgrid[[1]],xgrid[[2]])
       } else {
       bin[as.integer(levels(factor(ind[,1]))),as.integer(levels(factor(ind[,2]))),
       as.integer(levels(factor(ind[,3])))] <- table(ind[,1],ind[,2],ind[,3])
       xgrid <- list(seq(ry[1,1]+dry[1]/2/ngrid[1],ry[2,1]+dry[1]/2/ngrid[1],length=ngrid[1]),
                      seq(ry[1,2]+dry[2]/2/ngrid[2],ry[2,2]+dry[2]/2/ngrid[2],length=ngrid[2]),
                      seq(ry[1,3]+dry[3]/2/ngrid[3],ry[2,3]+dry[3]/2/ngrid[3],length=ngrid[3]))
       if(!is.null(fu)) u <- fu(xgrid[[1]],xgrid[[2]],xgrid[[3]])
       }
   }
   if(!is.null(fu)) u <- u*n*prod(dry/ngrid)
   dens <- laws(bin,qlambda=qlambda,model="Poisson",eta=eta,lkern=lkern,
                hinit=hinit,hincr=hincr,hmax=hmax,graph=graph,demo=demo,
                u=u,symmetric=symmetric)$theta
   dens <- dens/sum(dens)/prod(dry/ngrid)
list(bin=bin,dens=dens,xgrid=xgrid,call=args)
}
