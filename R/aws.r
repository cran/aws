############################################################################
#
# univariate Adaptive Weights Smoothing
#
# Copyright Weierstrass Instiute for Applied Analysis and Stochastics 
#           J. Polzehl 2000
############################################################################

awsuni <- function(y, lambda=3, gamma=1.3, eta =4, s2hat = NULL, kstar =
length(radii),radii = c(1:8,(5:12)*2,(7:12)*4,(7:12)*8,(7:10)*16,(6:8)*32,
         (5:8)*64,(5:8)*128,(5:8)*256),rmax=max(radii),
         graph = FALSE,z0 = NULL, eps = 1e-08, control="dyadic",demomode=FALSE)
{
# requires  dyn.load("aws.so") 
#
#   y - observed values (ordered by value of independent variable)
#   lambda - main smoothing parameter (should be approximately 3)
#   gamma  - allow for increase of variances (over minsk) by factor gamma
#   eta   - main control parameter (should be approximately 4)   
#   s2hat - initial variance estimate (if available,
#           can be either a number (homogeneous case), a vector of same length 
#           as y (inhomogeneous variance) or NULL (a homogeneous variance estimate
#           will be generated in this case)
#   kstar - number of iterations to perform (set to min(kstar, length(radii)))
#   radii - radii of neighbourhoods used
#   graph - logical, if TRUE progress (for each iteration) is illustrated grahically,
#           if FALSE the program runs until the final estimate is obtained 
#           (much faster !!!)
#   z0    - allows for submission of "true" values for illustration puposes only
#           if graph=TRUE  MSE and MAE are reported for each iteration step
#   eps - stop iteration if ||(yhatnew - yhat)||^2 < eps * sum(s2hat)
#   control - the control step is performed in either a dyadic sceme
#           ("dyadic") or using all previous estimates (otherwise)
#   demomode - only active if graph=TRUE, causes the program to wait after displaying the 
#           results of an iteration step
#
        radii <- radii[radii<=rmax]
        kstar <- min(kstar,length(radii))
        args <- list(lambda=lambda,gamma=gamma,eta=eta,s2hat = s2hat, 
                     kstar = kstar, radii=radii, rmax=rmax)
        if(graph) oldpar <- par(mfrow=c(1,3))
        ind <- trunc(radii)
        ind <- ind[ind>0]
        if(is.null(ind)) ind <- 1:kstar
        kstar <- min(kstar, length(ind))
        newcontr <- numeric(kstar)
        if(control=="dyadic") newcontr[2^(0:(log(kstar)/log(2)))] <- 1
        else newcontr[1:kstar] <- 1
        cat("Control sceme: ",newcontr,"\n")
        n <- length(y)
        x <- 1:n
        lam0 <- lambda
        lambda <- lambda^2
        gamma <- gamma^2
# generate a variance estimate if needed
        if(is.null(s2hat))
                s2hat <- (IQR(diff(y))/1.908)^2
# expand variance estimate in case of homogeneous variance
        if(length(s2hat) == 1)
                s2hat <- rep(s2hat, n)
# now initialize 
        yhat <- y
        sk <- skmin <- s2hat
        kern <- exp( - seq(0, 6, 0.3))
        lambda.3 <- lambda * 0.3
        controls <- numeric(2*n)
        dim(controls) <- c(2,n)
        controls[1,] <- y-eta*sqrt(s2hat)
        controls[2,] <- y+eta*sqrt(s2hat)
        if(ind[1] > 1) {
# nontrivial first neighbourhood (should only be used for small signal/noise)
                z <- .Fortran("locuini2",
                        as.integer(n),
                        as.single(y),
                        yhat = as.single(y),
                        sk = as.single(sk),
                        as.integer(ind[1] - 1),
                        as.single(numeric(2 * ind[1] - 1)),
                        as.single(numeric(2 * ind[1] - 1)),
                        as.single(s2hat),PACKAGE="aws")
                yhat <- z$yhat
                sk <- skmin <- z$sk
                if(graph) {
                        plot(x, y)
                        if(!is.null(z0)) lines(x, z0, col = 3)
                        lines(x, yhat, col = 2)
                        lines(x, yhat-sqrt(lambda*sk),col=4,lty=2)
                        lines(x, yhat+sqrt(lambda*sk),col=4,lty=2)
                        if(!is.null(z0)) lines(x, yhat, col = 2)
                        title(paste("Estimate  Iteration ", 0,
                        "  N(U) = ", 2 * ind[1] - 1))
                        ylim <- range(y - yhat)
                        if(!is.null(z0)) ylim <- range(ylim, z0 - yhat)
                        plot(x, y - yhat, ylim = ylim)
                        if(!is.null(z0)) lines(x, z0 - yhat, col = 3)
                        lines(x, yhat - yhat, col = 2)
                        title("Residuals")
                        plot(x, sqrt(sk))
                        title(paste("sigmahat\n l=", sqrt(lambda), "g=", sqrt(
                                gamma), "shat=", signif(sqrt(mean(s2hat)), 3)))
            if(demomode) {
            cat("press ENTER to continue")
            readline()
            }
                }
        }
        if(graph) {
           if(!is.null(z0))
                cat("Iteration ", 0, "MSE:", mean((yhat - z0)^2), "MAE:", 
                        mean(abs(yhat - z0)), "\n")
           for(k in 2:kstar) {
                yhatold <- yhat
                z <- .Fortran("locuniw",
                        as.integer(n),
                        as.single(y),
                        as.single(yhat),
                        yhat = as.single(yhat),
                        as.single(sk),
                        as.single(s2hat),
                        sk = as.single(sk),
                        controls=as.single(controls),
                        as.integer(newcontr[k]),
                        as.single(skmin),
                        as.integer(ind[k] - 1),
                        as.single(lambda.3),
                        as.single(eta),
                        as.single(gamma),
                        as.single(kern),PACKAGE="aws")[c("yhat","sk","controls")]
                yhat <- z$yhat
                sk <- z$sk
                controls <- z$controls
                skmin <- pmin(skmin, sk)
                if(!is.null(z0))
                cat("Iteration ", k-1, "MSE:", mean((yhat - z0)^2), "MAE:",
                        mean(abs(yhat - z0)), "\n")
                plot(x, y)
                if(!is.null(z0)) lines(x, z0, col = 3)
                lines(x, yhat, col = 2)
                lines(x, yhat-sqrt(lambda*sk),col=4,lty=2)
                lines(x, yhat+sqrt(lambda*sk),col=4,lty=2)
                if(!is.null(z0)) lines(x, yhat, col = 2)
                title(paste("Estimate  Iteration ", k-1, "  N(U) = ", 2 * ind[k] - 1))
                ylim <- range(y - yhat)
                if(!is.null(z0)) ylim <- range(ylim, z0 - yhat)
                plot(x, y - yhat, ylim = ylim)
                if(!is.null(z0)) lines(x, z0 - yhat, col = 3)
                lines(x, yhat - yhat, col = 2)
                title("Residuals")
                plot(x, sqrt(sk))
                title(paste("sigmahat    mean(shat)=", 
                      signif(sqrt(mean(sk)), 3)))
                if(sum((yhatold - yhat)^2) <= sum(eps * s2hat)) break
            if(demomode) {
            cat("press ENTER to continue")
            readline()
            }
                }
             par(oldpar)
        }
        else {
                z <- .Fortran("locunial",
                        as.integer(kstar),
                        as.integer(n),
                        as.single(y),
                        as.single(yhat),
                        yhat = as.single(yhat),
                        as.single(sk),
                        as.single(s2hat),
                        sk = as.single(sk),
                        controls=as.single(controls),
                        as.integer(newcontr),
                        as.single(skmin),
                        as.integer(ind),
                        as.single(lambda.3),
                        as.single(eta),
                        as.single(gamma),
                        as.single(kern),
                        as.single(sum(eps * s2hat)),PACKAGE="aws")[c("yhat","sk")]
                yhat <- z$yhat
                sk <- z$sk
        }
        if(!is.null(z0))
                cat("kstar=", kstar, "MSE:", mean((yhat - z0)^2), "MAE:",
                        mean(abs(yhat - z0)), "\n")
        list(yhat = yhat, shat = sqrt(sk),args=args)
}

############################################################################
#
# bivariate local constant smoothing
#
# Copyright Weierstrass Instiute for Applied Analysis and Stochastics
#           J. Polzehl 2000
############################################################################

awsbi <- function(y, lambda=3, gamma=1.3, eta = 4, 
     s2hat = NULL, kstar = length(radii), rmax=max(radii),
     radii=c((1:8)/2,4.4,5.,(6:10),(6:10)*2), graph = FALSE, 
     u0 = NULL,control="dyadic",demomode=FALSE, colors=gray((0:255)/255))
{
# requires  dyn.load("aws.so") 
#
#   y - observed values 
#   lambda - main smoothing parameter (should be approximately 3)
#   gamma  - allow for increase of variances (over minsk) by factor gamma
#   eta   - main control parameter (should be approximately 4)   
#   s2hat - initial variance estimate (if available,
#           can be either a number (homogeneous case), a matrix of same dimension  
#           as y (inhomogeneous variance) or NULL (a homogeneous variance estimate
#           will be generated in this case)
#   kstar - number of iterations to perform (set to min(kstar, length(radii)))
#   radii - radii of neighbourhoods used
#   graph - logical, if TRUE progress (for each iteration) is illustrated grahically,
#           if FALSE the program runs until the final estimate is obtained 
#           (much faster !!!)
#   colors - color sceme to be used for images
#   u0    - allows for submission of "true" values for illustration puposes only
#           if graph=TRUE  MSE and MAE are reported for each iteration step
#   control - the control step is performed in either a dyadic sceme
#           ("dyadic") or using all previous estimates (otherwise)
#   demomode - only active if graph=TRUE, causes the program to wait after displaying the 
#           results of an iteration step
#
        radii <- radii[radii<=rmax]
        kstar <- min(kstar,length(radii))
        args <- list(lambda=lambda,gamma=gamma,eta=eta,s2hat = s2hat, 
                     kstar = kstar, radii=radii)
        l2 <- r2 <- single(kstar)
        newcontr <- numeric(kstar)
        if(control=="dyadic") newcontr[2^(0:(log(kstar)/log(2)))] <- 1
        else newcontr[1:kstar] <- 1
        cat("Control sceme: ",newcontr,"\n")
        dy <- dim(y)
        if(is.null(dy)||length(dy)>2) stop("y should have dimension 2")
        nx <- dy[1]
        ny <- dy[2]
        n <- nx * ny
        if(graph) oldpar <- par(mfrow = c(1, 3))
        kiii <- radii^2
# get number of points in neighbourhoods
        iii <- getnubi(radii^2, c(1, 1))
        lambda <- lambda^2
        gamma <- gamma^2
        kern <- exp( - seq(0, 6, 0.3))
        lambda <- lambda * 0.3
# generate a variance estimate if needed
        if(length(s2hat) == 0) {
           s2hat <- (IQR(diff(y))/1.908)^2
           cat("Estimated variance:",s2hat,"\n")
           }
# expand s2hat in case of homogeneous variance
        if(length(as.vector(s2hat)) == 1) s2hat <- matrix(rep(s2hat, 
                        n), ncol = ny) 
        controls <- numeric(2*n)
        dim(controls) <- c(2,nx,ny)
        controls[1,,] <- y-eta*sqrt(s2hat)
        controls[2,,] <- y+eta*sqrt(s2hat)
        yhat <- y
        minsk <- sk <- s2hat
        if(!is.null(u0)) {
                l2[1] <- mean((yhat - u0)^2)
                r2[1] <- mean(abs(yhat - u0))
                cat("Iteration", 0, "nu=", iii[1], "MSE", l2[1], 
                                "MAE", r2[1], "\n")
        }
        kiiinit <- kiii[1]
        if(kiiinit > 5) kiiinit <- 5
# nontrivial first neighbourhood (should only be used for small signal/noise)
        z <- .Fortran("locbinis",
                as.integer(nx),
                as.integer(ny),
                as.single(y),
                yhat = as.single(yhat),
                as.single(s2hat),
                sk = as.single(sk),
                as.single(kiiinit + 0.001),PACKAGE="aws")
        yhat <- z$yhat
        sk <- z$sk
        if(graph) {
        for(k in 2:kstar) {
                z <- .Fortran("locbiw",
                        as.integer(nx),
                        as.integer(ny),
                        as.single(y),
                        as.single(yhat),
                        yhat = as.single(yhat),
                        as.single(sk),
                        sk = as.single(sk),
                        controls=as.single(controls),
                        as.integer(newcontr[k]),
                        as.single(minsk),
                        as.single(kiii[k] + 0.0001),
                        as.single(s2hat),
                        as.single(lambda),
                        as.single(eta),
                        as.single(gamma),
                        as.single(kern),PACKAGE="aws")[c("yhat","sk","controls")]
                yhat <- z$yhat
                sk <- z$sk
                controls <- z$controls
                minsk <- pmin(sk, minsk)
                if(!is.null(u0)) {
                        l2[k] <- mean((yhat - u0)^2)
                        r2[k] <- mean(abs(yhat - u0))
                }
                        image(matrix(y, ncol = ny),col=colors)
                        title("original image")
                        image(matrix(yhat, ncol = ny),zlim=range(y),col=colors)
                        title(paste("Estimate  Iteration ", k-1, "  N(U) = ", iii[k]))
                        image(matrix(log(sk), ncol = ny),col=colors)
                        title(paste("log(var(yhat))"," Mean Var:",signif(mean(sk),3)))
                if(!is.null(u0))
                        cat("Iteration", k-1, "nu=", iii[k], "MSE", l2[k],
                                "MAE", r2[k], "\n")
            if(demomode) {
            cat("press ENTER to continue")
            readline()
            }
                gc()
                }
                par(oldpar)
        list(yhat = matrix(yhat, ncol = ny), shat = matrix(sk, ncol = ny), 
             nu = iii, l2 = l2, r2 = r2, args=args)
        }
        else {
           z <- .Fortran("locbiall",
                        as.integer(kstar),
                        as.integer(nx),
                        as.integer(ny),
                        as.single(y),
                        as.single(yhat),
                        yhat = as.single(yhat),
                        as.single(sk),
                        sk = as.single(sk),
                        as.single(controls),
                        as.integer(newcontr),
                        as.single(minsk),
                        as.single(kiii),
                        as.single(s2hat),
                        as.single(lambda),
                        as.single(eta),
                        as.single(gamma),
                        as.single(kern),PACKAGE="aws")[c("yhat","sk")]
                if(!is.null(u0)){
                        l2[kstar] <- mean((z$yhat - u0)^2)
                        r2[kstar] <- mean(abs(z$yhat - u0))
                        cat("Iteration ", kstar-1, "nu=", iii[kstar], "MSE", l2[kstar],
                                "MAE", r2[kstar], "\n")
                        }
        list(yhat = matrix(z$yhat, ncol = ny),
             shat = matrix(z$sk, ncol = ny),
             nu = iii,  args=args)
        }
}

############################################################################
#
# trivariate local constant smoothing
#
# Copyright Weierstrass Instiute for Applied Analysis and Stochastics
#           J. Polzehl 2000
############################################################################


awstri <- function(y, lambda = 3, gamma = 1.3 , eta = 4, s2hat = NULL, 
    kstar = length(radii), rmax=max(radii), weight = c(1,1,1), 
    radii = c((1:4)/2,2.3,(5:12)/2,7:9,10.5,12,13.5),control="dyadic")
{
# requires  dyn.load("aws.so") 
#
#   y - observed values (ordered by value of independent variable)
#   lambda - main smoothing parameter (should be approximately 3)
#   gamma  - allow for increase of variances (over minsk) by factor gamma
#   eta   - main control parameter (should be approximately 4)   
#   s2hat - initial variance estimate (if available,
#           can be either a number (homogeneous case), a vector of same length 
#           as y (inhomogeneous variance) or NULL (a homogeneous variance estimate
#           will be generated in this case)
#   kstar - number of iterations to perform (set to min(kstar, length(radii)))
#   weight - excentricities of ellipsoids used as neighbourhoods 
#            used to weight distances in coordinate directions
#   radii - radii of neighbourhoods used
#   control - the control step is performed in either a dyadic sceme
#           ("dyadic") or using all previous estimates (otherwise)
#
# Speicherschonende Variante
# mit Fortran requires  dyn.load.shared("./image3.so")
        radii <- radii[radii<=rmax]
        kstar <- min(kstar,length(radii))
        args <- list(lambda=lambda,gamma=gamma,eta=eta,s2hat = s2hat, 
                     kstar = kstar, radii=radii)
        dy <- dim(y)
        if(length(dy) != 3) stop("y is not a 3-dimensional array")
        if(is.null(weight)) weight <- rep(1, 3)
        if(is.null(radii)) stop("No neigborhood defined")
        radii2 <- radii^2
        nx <- dy[1]
        ny <- dy[2]
        nz <- dy[3]
        newcontr <- numeric(kstar)
        if(control=="dyadic") newcontr[2^(0:(log(kstar)/log(2)))] <- 1
        else newcontr[1:kstar] <- 1
    cat("Control sceme: ", newcontr,"\n")
        lambda <- lambda^2
        gamma <- gamma^2
#         kern <- exp( - seq(0, 6, 0.3)/1.44)
        kern <- exp( - seq(0, 6, 0.3))
        lambda <- lambda*.3
        n <- nx * ny * nz
        if(length(s2hat) == 0) s2hat <- (IQR(diff(y))/1.908)^2
        if(length(as.vector(s2hat)) == 1) homogeneous <- TRUE
        #now precompute neighbourhoods
        yhat <- y
        if(homogeneous) minsk <- sk <- array(s2hat,dim(y))
    else minsk <- sk <- s2hat
        controls <- numeric(2*n)
        dim(controls) <- c(2,nx,ny,nz)
        controls[1,,,] <- y-eta*sqrt(s2hat)
        controls[2,,,] <- y+eta*sqrt(s2hat)
        if(homogeneous){
        z <- .Fortran("loctria0",
                as.integer(kstar),
                as.integer(nx),
                as.integer(ny),
                as.integer(nz),
                as.single(y),
                as.single(yhat),
                yhat = as.single(yhat),
                as.single(sk),
                sk = as.single(sk),
                as.single(controls),
                as.integer(newcontr),
                as.single(radii2),
                as.single(s2hat),
                as.single(lambda),
                as.single(eta),
                as.single(weight),
                as.single(kern),PACKAGE="aws")[c("yhat","sk")]
                }
        else{
        z <- .Fortran("loctrial",
                      as.integer(kstar),
                      as.integer(nx),
                      as.integer(ny),
                      as.integer(nz),
                      as.single(y),
                      as.single(yhat),
                      yhat = as.single(yhat),
                      as.single(sk),
                      sk = as.single(sk),
                      as.single(controls),
                      as.integer(newcontr),
                      as.single(minsk),
                      as.single(radii2),
                      as.single(s2hat),
                      as.single(lambda),
                      as.single(eta),
                      as.single(gamma),
                      as.single(weight),
                      as.single(kern),PACKAGE="aws")[c("yhat","sk")]
                      }
    list(yhat = array(z$yhat, dim = dy), shat = array(z$sk,
         dim = dy), args = args)
}


############################################################################
#
# get number of pixels in bivariate neighbourhoods
#
############################################################################

getnubi <- function(radiusq, weights)
{
        nu <- numeric(length(radiusq))
        nu <- .Fortran("getnubi",
                as.single(radiusq),
                as.single(weights),
                nu = as.integer(nu),
                as.integer(length(radiusq)),PACKAGE="aws")$nu
        nu
}

############################################################################
#
# get number of pixels in trivariate neighbourhoods
#
############################################################################

getnutri <- function(radiusq, weights)
{
        nu <- nu3 <- numeric(length(radiusq))
        nu3 <- .Fortran("getnubi",
                as.single(radiusq),
                as.single(weights[1:2]),
                nu = as.integer(nu3),
                as.integer(length(radiusq)),PACKAGE="aws")$nu
        for( i in 1:trunc(max(sqrt(radiusq))/weights[3])) {
           radius2 <- radiusq-i^2*weights[3]
           ind <- (1:length(radiusq))[radius2>=0]
           if(length(ind)<1) break
           nu3[ind] <- nu3[ind] + 2*.Fortran("getnubi",
                         as.single(radius2[ind]),
                         as.single(weights[1:2]),
                         nu = as.integer(nu[ind]),
                         as.integer(length(ind)),PACKAGE="aws")$nu
        }
        nu3
}


.First.lib <- function(lib, pkg) {
  if(version$major==0)
    stop("This version for R 1.00 or later")
  library.dynam("aws", pkg, lib)
}
