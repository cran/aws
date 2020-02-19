### R code from vignette source 'aws-Example.Rnw'

###################################################
### code chunk number 1: 0
###################################################
options(digits=3)


###################################################
### code chunk number 2: 1a
###################################################
library(aws)
fofx1 <- c(rep(0,25),rep(-1,20),rep(1,20), rep(-2,10),rep(2,5), 
           rep(-1,25),rep(-.5,30),rep(0,35)) 
set.seed(1)
y1 <- rnorm(fofx1,fofx1,.3)


###################################################
### code chunk number 3: 1b
###################################################
u1 <- matrix(0,64,64)
ind0 <- seq(0,1,length=64)
ind <- outer(ind0^2,ind0^2,"+")
u1[ind  > .95] <- u1[ind >.95] + 2
u1[ind < .6] <- u1[ind < .6] -2
u1[ind < .35] <- u1[ind < .35] +3
u1[ind < .15] <- u1[ind < .15] -2
u1[ind < .05] <- u1[ind < .05] +3
u1 <- u1*(1-2*outer(ind0,ind0,">"))
z1 <- u1+rnorm(u1)


###################################################
### code chunk number 4: 1c
###################################################
u2 <- u1+5*ind
z2 <- u2+rnorm(u1)


###################################################
### code chunk number 5: 2a
###################################################
yhat0 <- kernsm(y1, h=10)


###################################################
### code chunk number 6: 2a
###################################################
yhat1 <- aws(y1, hmax=100)


###################################################
### code chunk number 7: 2b
###################################################
par(mfrow=c(1,3), mar=c(3,3,3,1), mgp=c(2,1,0))
plot(y1)
lines(yhat1@theta, col=2)
lines(fofx1, col=3)
title("AWS estimate")
plot(yhat1@ni)
title("Sum of weights")
plot(y1)
lines(kernsm(y1,.609)@yhat, col=2)
lines(fofx1, col=3)
title("MSE optimal kernel estimate")


###################################################
### code chunk number 8: 2c
###################################################
setCores(2)
zhat1a <- aws(z1, hmax=8)
zhat1b <- paws(z1, hmax=10, patchsize=1)


###################################################
### code chunk number 9: 2d
###################################################
par(mfrow=c(2,3), mar=c(3,3,3,1), mgp=c(2,1,0))
image(z1, col=grey(0:255/255))
title("Noisy original")
image(zhat1a@theta, col=grey(0:255/255))
title("AWS reconstruction")
image(zhat1a@ni, col=grey(0:255/255))
title("AWS sum of weights")
image(u1, col=grey(0:255/255))
title("True image")
image(zhat1b@theta, col=grey(0:255/255))
title("PAWS reconstruction")
image(zhat1b@ni, col=grey(0:255/255))
title("PAWS sum of weights")


###################################################
### code chunk number 10: 2e
###################################################
zhat2a <- aws(z2, hmax=8)
zhat2b <- paws(z2, hmax=10)


###################################################
### code chunk number 11: 2f
###################################################
par(mfrow=c(2,3), mar=c(3,3,3,1), mgp=c(2,1,0))
image(z2, col=grey(0:255/255))
title("Noisy original")
image(zhat2a@theta, col=grey(0:255/255))
title("AWS reconstruction")
image(zhat2a@ni, col=grey(0:255/255))
title("AWS sum of weights")
image(u2, col=grey(0:255/255))
title("True image")
image(zhat2b@theta, col=grey(0:255/255))
title("PAWS reconstruction")
image(zhat2b@ni, col=grey(0:255/255))
title("PAWS sum of weights")


###################################################
### code chunk number 12: 3a
###################################################
zhat1c <- kernsm(z1,.9)@yhat
zhat1d <- ICIsmooth(z1, hmax=8, thresh=.8, presmooth=TRUE)@yhat
zhat1e <- ICIcombined(z1, hmax=8, nsector=8, thresh=.8, 
                      presmooth=TRUE)@yhat


###################################################
### code chunk number 13: 3b
###################################################
par(mfrow=c(1,4), mar=c(3,3,3,1), mgp=c(2,1,0))
image(z1, col=grey(0:255/255))
title("Noisy original")
image(zhat1c, col=grey(0:255/255))
title("optimal kernel estimate")
image(zhat1d, col=grey(0:255/255))
title("adaptation over h")
image(zhat1e, col=grey(0:255/255))
title("adaptation over h and sectorial")


###################################################
### code chunk number 14: 3a
###################################################
zhat1f <- nlmeans(z1, .85, 1, searchhw=6)$theta 


###################################################
### code chunk number 15: 4a
###################################################
zhat1f <- TV_denoising(z1, .93)
zhat1g <- TGV_denoising(z1, .92, 4)


###################################################
### code chunk number 16: 3b
###################################################
par(mfrow=c(1,4), mar=c(3,3,3,1), mgp=c(2,1,0))
image(z1, col=grey(0:255/255))
title("Noisy original")
image(zhat1e, col=grey(0:255/255))
title("NL-Means estimate")
image(zhat1f, col=grey(0:255/255))
title("Optimal TV reconstruction")
image(zhat1g, col=grey(0:255/255))
title("Optimal TGV reconstruction")


