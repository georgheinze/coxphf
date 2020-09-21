coxphfplot <-
  function(
    formula, 
    data, 
    profile,      # righthand formula
    pitch=.05,    # distances between points in std's
    limits,       # vector of MIN & MAX in std's, default=extremes of both CI's +-0.5 std. of beta
    alpha=.05, 
    maxit = 50, 
    maxhs = 5, 
    epsilon = 1e-006, 
    maxstep = 0.5, 
    firth = TRUE, 
    penalty=0.5,
    adapt=NULL,
    legend="center",	# see R-help of function <legend>, "" produces no legend
    ...			# other parameters to <legend>
  ) {
    ### by MP und GH, 2006-10
    ###
    
    call <- match.call()
    mf <- match.call(expand.dots =FALSE)
    m <- match(c("formula","data"), names(mf), 0L)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    ## call coxphf:
    fit <- coxphf(formula=formula, data=data, alpha=alpha, maxit=maxit, maxhs=maxhs, 
                  epsilon=epsilon, maxstep=maxstep, firth=firth, pl=TRUE, penalty=penalty, adapt=adapt)
    coefs <- coef(fit)           
    covs <- fit$var              
    n <- nrow(data)
    
    obj <- decomposeSurv(formula, data, sort=TRUE)
    NTDE <- obj$NTDE
    mmm <- cbind(obj$mm1, obj$timedata)
    
    cov.name <- obj$covnames
    
    k <- ncol(obj$mm1)          # number covariates
    ones <- matrix(1, n, k+NTDE)
    
    ## standardise
    sd1 <- apply(as.matrix(obj$mm1),2,sd)
    sd2 <- apply(as.matrix(obj$timedata),2,sd)
    Z.sd <- c(sd1, sd2 * sd1[obj$timeind])
    obj$mm1 <- scale(obj$mm1, FALSE, sd1)
    obj$timedata <- scale(obj$timedata, FALSE, sd2)
    mmm <- cbind(obj$mm1, obj$timedata)
    
    CARDS <- cbind(obj$mm1, obj$resp, ones, obj$timedata)   
    PARMS <- c(n, k, firth, maxit, maxhs, maxstep, epsilon, 1, 0.0001, 0, 0, 0, 0, NTDE, penalty)
    
    #--> nun Berechnungen fuer Schleife
    formula2 <- as.formula(paste(as.character(formula)[2], 
                                 as.character(profile)[2], sep="~"))
    obj2 <- decomposeSurv(formula2, data, sort = TRUE)
    cov.name2 <- obj2$covnames       # Labels der Test-Fakt.
    pos <- match(cov.name2, cov.name) ## Position der Testfakt.
    
    std.pos <- diag(fit$var)[pos]^.5
    if(missing(limits)) {
      lim.pl <- (c(log(fit$ci.lower[pos]), log(fit$ci.upper[pos])) - 
                   coef(fit)[pos]) / std.pos
      limits <- c(min(qnorm(alpha / 2), lim.pl[1]) - .5,
                  max(qnorm(1 - alpha / 2), lim.pl[2]) + .5)
    }
    limits <- c(floor(limits[1]/pitch)*pitch, ceiling(limits[2]/pitch)*pitch)
    knots <- seq(limits[1], limits[2], pitch)
    iflag <- rep(1, k+NTDE)
    if(!is.null(adapt)) iflag<-adapt
    iflag[pos] <- 0
    offset <- rep(0, k+NTDE)
    nn <- length(knots)
    res <- matrix(knots, nn, 3) #initialisiere Werte
    dimnames(res) <- list(1:nn, c("std", cov.name2, "log-likelihood"))
    for(i in 1:nn) {
      res[i, 2] <- coefs[pos] + covs[pos,pos]^.5 * knots[i]
      offset[pos] <- res[i, 2] * Z.sd[pos]  
      IOARRAY <- rbind(iflag, offset, matrix(0, 1+k+NTDE, k + NTDE))
      
      # --------------- Aufruf Fortran - Makro FIRTHCOX ----------------------------------
      value <- .Fortran("firthcox",
                        CARDS,
                        outpar = PARMS,
                        outtab = IOARRAY, PACKAGE="coxphf")
      if(value$outpar[8])
        warning("Error in routine FIRTHCOX; parms8 <> 0")
      res[i, 3] <- value$outpar[11]
    }
    
    #### Graphischer Output:
    my.par <- act.par <- par()
    my.par$mai[3] <- 1.65 * act.par$mai[3]
    ###if(legend != "") 
    ###	my.par$mai[1] <- 2.00 * act.par$mai[1]
    par(mai=my.par$mai)
    ind <- (1:nn)[round(4*res[,1])==round(4*res[,1], 10)]
    if (length(ind)==0) ind <- 1:nn
    pp <- max(res[,3]) - .5*res[,1]^2
    plot(res[,-1], type="l", xlab=paste("BETA of", cov.name2)) ##Profile likelihood
    ##lines(res[,2], pp, lty=4)  #<<<Wald approximative profile lik. >>>
    points(res[res[,1]==0,2], max(res[,3])) ##Maximum of likelihood
    segments(min(res[,2]), max(res[,3])-.5*qchisq(1-alpha, 1), 
             max(res[,2]), max(res[,3])-.5*qchisq(1-alpha, 1), lty=3) ##refer.line
    yy <- par("usr")[4] - (par("usr")[4]-par("usr")[3])*c(.9, .95)
    segments(coef(fit)[pos] - qnorm(alpha/2) * std.pos, yy[1], 
             coef(fit)[pos] - qnorm(1 - alpha/2) * std.pos, yy[1], lty=6) ##Wald-CI
    segments(log(fit$ci.lower[pos]), yy[2], 
             log(fit$ci.upper[pos]), yy[2], lty=8) ## prof.pen.lik.-CI
    axis(side=3, at=res[ind,2], labels=res[ind,1])
    mtext("distance from maximum in standard deviations", side=3, line=3)
    if(legend != "")
      legend(legend,
             ##x=par("usr")[1], 
             ##y=par("usr")[3]-.24*(par("usr")[4] - par("usr")[3]), 
             legend=c("Profile penalized likelihood", 
                      paste(100*(1-alpha), "% - reference line", sep=""), 
                      "Maximum of Likelihood",
                      "Wald - confidence interval", 
                      "Profile penalized likelihood CI"),
             lty=c(1, 3, -1, 6, 8),
             pch=c(-1, -1, 1, -1, -1),
             bty="n", ...
      )
    par(mai=act.par$mai)
    title("Profile likelihood")
    invisible(res)
  }