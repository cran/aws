fw3<-1:10
n3<-10
as.vector(outer(c(1,.8,.6,.4,.2),fw3[-n3],"*")+outer(1-c(1,.8,.6,.4,.2),fw3[-1],"*"))


library(aws)
y<-rnorm(u,u)
#debug(hiraraws)
yhat<-hiraraws(y,hmax=100,qlambda=.995,graph=TRUE)


library(aws)
x<-seq(0,2*pi,.001)
y<-sin(x)+rnorm(x,0,.1)
yhat<-hiraraws(y,hmax=1000,qlambda=.995,graph=TRUE)

library(aws)
x<-seq(0,1,.001)
y <- 1000*x^3+100*x^2+10*x+rnorm(x)
yhat<-hiraraws(y,hmax=1000,qlambda=.995,graph=TRUE)

library(aws)
y <- rnorm(10000)
yhat<-hiraraws(y,hmax=1000,qlambda=.995,graph=TRUE)

library(aws)
y <- array(rnorm(100^3),c(100,100,100))
yhat <- aws(y,hmax=10,lseq=c(1.9,1.5,1.3,1.3,1.3,1.3,rep(1.1,8)),u=0,testprop=TRUE)

alpha=.1: lseq=c(1.9,1.45,1.15,1.05,1.25,1.3,1.05,1,.95,.95,1,1.1,1.15,rep(1,20))
alpha=.2:.75*c(1.4,1.25,1.05,1.0,1.15,1.25,1.0,1,.95,.95,.95,1.05,1.1,.95,rep(1,20))

c(1.9,1.45,1.15,1.05,1.25,1.3,1.05,1,.95,.95,1,1.1,1.15,rep(1,20))[1:30]/(.75*c(1.4,1.25,1.05,1.0,1.15,1.25,1.0,1,.95,.95,.95,1.05,1.1,.95,rep(1,20)))[1:30]


Parameter Sequenzen fuer Gausssche Modelle: 

3D:  alpha=0.2,  plateu=0.25,  lambda=11.1
 
Korrekturfaktoren fuer lambda: lseq <- c(1.,1.,0.83,0.83,1.,1.05,0.86,0.86,0.86,0.86,0.86,1.,1.,0.86,0.9,0.95,1,1,...)

2D:  alpha=0.15,  plateu=0.25,  lambda=12.7

Korrekturfaktoren fuer lambda: lseq <- c(1,.88,0.75,0.95,.8,.75,0.85,1,0.86,0.86,0.95,.925,.9,0.925,0.925,0.95,1,1,...)

lseq[k] ist der Korrekturfaktor bei Bandweite hakt=1.25^(k/d)  

#aws_1.4-2.tar.gz 
library(aws)
set.seed(1)
y <- array(rnorm(512^2),c(512,512))
lseq <- c(1,.88,0.75,0.95,.8,.75,0.85,1,0.86,0.86,0.95,.925,.9,0.925,0.925,0.95,1,1)
yhat <- aws(y,hmax=10,lseq=lseq,u=0,testprop=TRUE,spmin=.25,qlambda=.94)

y <- array(rnorm(50^3),c(50,50,50))
lseq <- c(1.,1.,0.83,0.83,1.,1.05,0.86,0.86,0.86,0.86,0.86,1.,1.,0.86,0.9,0.95,1,1)
yhat <- aws(y,hmax=10,lseq=lseq,u=0,testprop=TRUE,spmin=.25,qlambda=.92)

