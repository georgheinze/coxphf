#' Plot the Penalized Profile Likelhood Function
#' 
#' Plots the penalized profile likelihood for a specified parameter.
#' 
#' This function plots the profile (penalized) log likelihood of the specified parameter. A symmetric shape of 
#' the profile (penalized) log likelihood (PPL) function allows use of Wald intervals, while an
#' asymmetric shape demands profile (penalized) likelihood intervals (Heinze & Schemper (2001)).
#'
#' @param formula a formula object, with the response on the left of the  operator, and the
#' model terms on the right. The response must be a survival object as returned by the 'Surv' function.
#' @param data a data.frame in which to interpret the variables named in the 'formula' argument.
#' @param profile a righthand formula specifying the plotted parameter, interaction or general term, e.g. \code{~ A} or \code{~ A : C}.
#' @param pitch distances between the interpolated points in standard errors of the parameter estimate, the default value is 0.05.
#' @param limits the range of the x-axis in terms of standard errors from the parameter estimate. The default values 
#' are the extremes of both confidence intervals, Wald and PL, plus or minus half a
#' standard error, respectively.
#' @param alpha the significance level (1-\eqn{\alpha} the confidence level, 0.05 as default).
#' @param maxit maximum number of iterations (default value is 50)
#' @param maxhs maximum number of step-halvings per iterations (default value is 5). 
#' The increments of the parameter vector in one Newton-Rhaphson iteration step are halved, 
#' unless the new likelihood is greater than the old one, maximally doing \code{maxhs} halvings.
#' @param epsilon specifies the maximum allowed change in penalized log likelihood to declare convergence. Default value is 0.0001.
#' @param maxstep specifies the maximum change of (standardized) parameter values allowed in one iteration. Default value is 2.5.
#' @param firth use of Firth's penalized maximum likelihood (\code{firth=TRUE}, default) or the standard maximum likelihood method (\code{firth=FALSE}) for fitting the Cox model.
#' @param penalty optional: specifies a vector of 1s and 0s, where 0 means that the corresponding parameter is fixed at 0, while 1 enables 
#' parameter estimation for that parameter. The length of adapt must be equal to the number of parameters to be estimated.
#' @param adapt strength of Firth-type penalty. Defaults to 0.5.
#' @param legend if FALSE, legends in the plot would be omitted (default is TRUE).
#' @param ... other parameters to legend
#'
#' @return A matrix of dimension \eqn{m \times 3}, with \eqn{m = 1/\code{pitch} + 1}. 
#' With the default settings, \eqn{m=101}. The column headers are:
#' \item{std}{the distance from the parameter estimate in standard errors}
#' \item{x}{the parameter value}
#' \item{log-likelihood}{the profile likelihood at \code{x}}
#' 
#' @export
#'
#' @examples
#' library(survival)
#' time<-c(1,2,3)
#' cens<-c(1,1,1)
#' x<-c(1,1,0)
#' sim<-cbind(time,cens,x)
#' sim<-data.frame(sim)
#' profplot<-coxphfplot(sim, formula=Surv(time,cens)~x, profile=~x)
#' 
#' @references Firth D (1993). Bias reduction of maximum likelihood estimates. \emph{Biometrika} 80:27--38.
#' 
#' Heinze G and Schemper M (2001). A Solution to the Problem of Monotone Likelihood in Cox Regression. \emph{Biometrics} 57(1):114--119. 
#' 
#' Heinze G (1999). Technical Report 10/1999: The application of Firth's procedure to Cox and logistic regression. Section of Clinical Biometrics, Department of Medical Computer Sciences, University of Vienna, Vienna.
#' @author Georg Heinze and Meinhard Ploner
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