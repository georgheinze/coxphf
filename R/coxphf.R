#' Cox Regression with Firth's Penalized Likelihood
#' 
#' Implements Firth's penalized maximum likelihood bias reduction method  for Cox regression 
#' which has been shown to provide a solution in case of monotone likelihood (nonconvergence of likelihood function). 
#' The program fits profile penalized likelihood confidence intervals which were proved to outperform 
#' Wald confidence intervals.
#' 
#' The phenomenon of monotone likelihood in a sample causes parameter estimates of a Cox model to diverge, with 
#' infinite standard errors. Therefore, classical maximum likelihood analysis fails; the usual Wald confidence 
#' intervals cover the whole range of real numbers. Monotone likelihood appears if there is single covariate
#' or a linear combination of covariates such that at each event time, out of all individuals being at risk at 
#' that time, the individual with the highest (or at each event time the individual with the lowest) value for that covariate or linear combination experiences the event. It was shown that 
#' analysis by Firth's penalized likelihood method, particularly in conjunction with the computation
#' of profile likelihood confidence intervals and penalized likelihood ratio tests is superior to maximum
#' likelihood analysis. It completely removes the convergence problem mentioned in the paragraph on CONVERGENCE of
#' the description of the function \code{coxph}. The \code{formula} may involve time-dependent effects or 
#' time-dependent covariates. The response may
#' be given in counting process style, but it cannot be used for multivariate failure times, as the program has no option
#' to fit a robust covariance matrix. The user is responsible for the independency of observations within each risk set, i.e., 
#' the same individual should not appear twice within the same risk set.
#' 
#' @param formula a formula object, with the response on the left and the model terms on the right. 
#' The response must be a survival object as returned by the 'Surv' function (see its documentation in the survival package)
#' @param data a data.frame in which to interpret the variables named in the 'formula' argument.
#' @param pl specifies if confidence intervals and tests should be based on the profile penalized log likelihood (\code{pl=TRUE}, the default) or on the Wald method (\code{pl=FALSE}).
#' @param alpha the significance level (1-\eqn{\alpha} = the confidence level), 0.05 as default.
#' @param maxit maximum number of iterations (default value is 50)
#' @param maxhs maximum number of step-halvings per iterations (default value is 5). 
#' The increments of the parameter vector in one Newton-Rhaphson iteration step are halved,
#' unless the new likelihood is greater than the old one, maximally doing \code{maxhs} halvings.
#' @param epsilon specifies the maximum allowed change in standardized parameter estimates to declare convergence. Default value is 1e-6.
#' @param gconv specifies the maximum allowed absolute value of first derivative of likelihood to declare convergence. Default value is 0.0001.
#' @param maxstep specifies the maximum change of (standardized) parameter values allowed in one iteration. Default value is 0.5.
#' @param firth use of Firth's penalized maximum likelihood (\code{firth=TRUE}, default) or the standard maximum likelihood method (\code{firth=FALSE}) for fitting the Cox model.
#' @param adapt optional: specifies a vector of 1s and 0s, where 0 means that the corresponding parameter is fixed at 0, while 1 enables parameter estimation for that parameter. The length of adapt must be equal to the number of parameters to be estimated.
#' @param penalty strength of Firth-type penalty. Defaults to 0.5.
#' 
#' @return The object returned is of the class \code{coxphf} and has the following attributes:
#'  \item{coefficients}{the parameter estimates}
#' \item{alpha}{the significance level = 1 - confidence level}
#' \item{var}{the estimated covariance matrix}
#' \item{df}{the degrees of freedom}
#' \item{loglik}{the null and maximimized (penalized) log likelihood}
#' \item{method.ties}{the ties handling method}
#' \item{iter}{the number of iterations needed to converge}
#' \item{n}{the number of observations}
#' \item{y}{the response}
#' \item{formula}{the model formula}
#' \item{means}{the means of the covariates}
#' \item{linear.predictors}{the linear predictors}
#' \item{method}{the estimation method (Standard ML or Penalized ML)}
#' \item{method.ci}{the confidence interval estimation method (Profile Likelihood or Wald)}
#' \item{ci.lower}{the lower confidence limits}
#' \item{ci.upper}{the upper confidence limits}
#' \item{prob}{the p-values}
#' \item{call}{the function call}
#' \item{iter.ci}{the numbers of iterations needed for profile likelihood confidence interval estimation, and for maximizing the restricted likelihood for p-value computation.}
#' 
#' @export
#'
#' @encoding UTF-8
#' @examples
#' # fixed covariate and monotone likelihood
#' library(survival)
#' time<-c(1,2,3)
#' cens<-c(1,1,1)
#' x<-c(1,1,0)
#' sim<-cbind(time,cens,x)
#' sim<-data.frame(sim)
#' coxphf(sim, formula=Surv(time,cens)~x) #convergence attained!
#' #coxph(sim, formula=Surv(time,cens)~x)  #no convergence!
#' # time-dependent covariate
#' test2 <- data.frame(list(start=c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
#'                          stop =c(2, 3, 6, 7, 8, 9, 9, 9,14,17),
#'                          event=c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
#'                          x    =c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0) ))
#' 
#' summary( coxphf( formula=Surv(start, stop, event) ~ x, pl=FALSE, data=test2))
#' 
#' 
#' # time-dependent effect
#' # the coxphf function can handle interactions of a (fixed or time-dependent)
#' # covariate with time
#' # such that the hazard ratio can be expressed as a function of time
#' 
#' summary(coxphf(formula=Surv(start, stop, event)~x+x:log(stop), data=test2, pl=FALSE, firth=TRUE))
#' 
#' # note that coxph would treat x:log(stop) as a fixed covariate
#' # (computed before the iteration process)
#' # coxphf treats x:log(stop) as a time-dependent covariate which changes (
#' # for the same individual!) over time
#' 
#' 
#' # time-dependent effect with monotone likelihood
#' 
#' test3 <- data.frame(list(start=c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
#'                          stop =c(2, 3, 6, 7, 8, 9, 9, 9,14,17),
#'                          event=c(1, 0, 0, 1, 0, 1, 1, 0, 0, 0),
#'                          x    =c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0) ))
#' 
#' summary( coxphf( formula=Surv(start, stop, event) ~ x+x:log(stop), pl=FALSE, maxit=400, data=test3))
#' 
#' 
#' # no convergence if option "firth" is turned off: 
#' # summary( coxphf(formula=Surv(start, stop, event) ~ x+x:log(stop), pl=F,
#' #                 data=test3, firth=FALSE)
#' 
#' 
#' data(breast)
#' fit.breast<-coxphf(data=breast, Surv(TIME,CENS)~T+N+G+CD)
#' summary(fit.breast)
#' 
#' @author Georg Heinze and Meinhard Ploner
#' @references Firth D (1993). Bias reduction of maximum likelihood estimates. \emph{Biometrika} 80:27--38.
#' 
#' Heinze G and Schemper M (2001). A Solution to the Problem of Monotone Likelihood in Cox Regression. \emph{Biometrics} 57(1):114--119. 
#' 
#' Heinze G (1999). Technical Report 10/1999: The application of Firth's procedure to Cox and logistic regression. Section of Clinical Biometrics, Department of Medical Computer Sciences, University of Vienna, Vienna.
#' @seealso [coxphfplot, coxphftest]
coxphf <-
function(
 formula,	# formula where .time. represents time for time-interactions
 data,
 pl=TRUE,
 alpha=0.05,
 maxit=50,
 maxhs=5,
 epsilon=1e-6,
 gconv=1e-4,
 maxstep=0.5,
 firth=TRUE,
 adapt=NULL,
 penalty=0.5
){
### by MP und GH, 2006-2018
  call <- match.call()
  mf <- match.call(expand.dots =FALSE)
  m <- match(c("formula","data"), names(mf), 0L)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  
  y <- model.response(mf)
  n <- nrow(data)
  x <- model.matrix(mt, mf)

  if(!is.logical(firth)) stop("Please set option firth to TRUE or FALSE.\n")
  if(!is.logical(pl)) stop("Please set option pl to TRUE or FALSE.\n")

	obj <- decomposeSurv(formula, data, sort=TRUE)
  NTDE <- obj$NTDE
	mmm <- cbind(obj$mm1, obj$timedata)
        
	cov.name <- obj$covnames
  k <- ncol(obj$mm1)          # number of covariates
  ones <- matrix(1, n, k+NTDE)
        
  if(!is.null(adapt)){
      if(k+NTDE != length(adapt)) stop("length of adapt must match the number of parameters to be estimated.")
      if(any(adapt !=0 & adapt !=1)) stop("adapt must consist of 0s and 1s exclusively.")
  }
          
	  ## standardise
	sd1 <- apply(as.matrix(obj$mm1),2,sd)
	sd2 <- apply(as.matrix(obj$timedata),2,sd)
	Z.sd <- c(sd1, sd2 * sd1[obj$timeind])
	obj$mm1 <- scale(obj$mm1, FALSE, sd1)
	obj$timedata <- scale(obj$timedata, FALSE, sd2)
	mmm <- cbind(obj$mm1, obj$timedata)
	   
	CARDS <- cbind(obj$mm1, obj$resp, ones, obj$timedata)	  
  PARMS <- c(n, k, firth, maxit, maxhs, maxstep, epsilon, 1, gconv, 0, 0, 0, 0, NTDE, penalty)
  IOARRAY <- rbind(rep(1, k+NTDE), matrix(0, 2+k+NTDE, k + NTDE))
  if(!is.null(adapt)) IOARRAY[1,]<-adapt
  if(NTDE>0)
  IOARRAY[4, (k+1):(k+NTDE)] <- obj$timeind
  storage.mode(CARDS) <- "double"
  storage.mode(PARMS) <- "double"
  storage.mode(IOARRAY) <- "double"

  ## --------------- Call Fortran routine FIRTHCOX ----------------------------------
  value <- .Fortran("firthcox",
                CARDS,
                outpar = PARMS,
                outtab = IOARRAY, PACKAGE="coxphf")
  if(value$outpar[8]) warning("Numerical problem in parameter estimation; check convergence.\n")
  outtab <- matrix(value$outtab, nrow=3+k+NTDE) #

  # --------------- create output object -------
  coef.orig <- outtab[3,  ]
  coefs <- coef.orig / Z.sd
  covs <- matrix(outtab[4:(k+3+NTDE), ], ncol=k+NTDE) / (Z.sd %*% t(Z.sd))

  dimnames(covs) <- list(cov.name, cov.name)
  vars <- diag(covs)
  names(coefs) <- cov.name
  df<-k+NTDE
  if(!is.null(adapt)) df<-sum(adapt>0)
  fit <- list(coefficients = coefs, alpha = alpha, var = covs, df = df,
              loglik = value$outpar[12:11], iter = value$outpar[10],
              method.ties = "breslow", n = n, ##terms = terms(formula),
              y = obj$resp, formula = formula, call = match.call())
  fit$means <- apply(mmm, 2, mean)
  fit$linear.predictors <- as.vector(scale(mmm, fit$means, scale=FALSE) %*% coefs)
	if(fit$iter>=maxit) warning("Convergence for parameter estimation not attained in ", maxit, " iterations.\n")
  if(firth) fit$method <- "Penalized ML"
  else fit$method <- "Standard ML"
  if(penalty != 0.5) fit$method<-paste(fit$method, " (penalty=",penalty,")",sep="")

  # --------------- Call Fortran routine PLCOMP ------------------------------------
  if(pl) {
    PARMS <- c(PARMS[1:7], qchisq(1-alpha, 1), gconv, 0, 0, 0, 0, NTDE, penalty)
    IOARRAY <- rbind(rep(1, k+NTDE), rep(0, k+NTDE), coef.orig, matrix(0, 6, k+NTDE))
    if(!is.null(adapt)) IOARRAY[1,]<-adapt
    if(NTDE>0) IOARRAY[4,(k+1):(k+NTDE)] <- obj$timeind
    storage.mode(PARMS) <- "double"
    storage.mode(IOARRAY) <- "double"
    value <- .Fortran("plcomp",
                        CARDS,
                        outpar = PARMS,
                        outtab = IOARRAY, PACKAGE="coxphf")
    if(value$outpar[9]) warning("Numerical problem in estimating confidence intervals; check convergence.\n")
    fit$method.ci <- "Profile Likelihood"
    fit$ci.lower <- exp(value$outtab[4,  ] / Z.sd)
    fit$ci.upper <- exp(value$outtab[5,  ] / Z.sd)
    fit$prob <- 1 - pchisq(value$outtab[6,  ], 1)
    fit$iter.ci<-t(value$outtab[7:9,])
    colnames(fit$iter.ci)<-c("Lower", "Upper", "P-value")
    rownames(fit$iter.ci)<-cov.name
    fit$ci.lower[fit$iter.ci[,1]>=maxit]<-NA
    fit$ci.upper[fit$iter.ci[,2]>=maxit]<-NA
    fit$prob[fit$iter.ci[,3]>=maxit]<-NA
    if(any(fit$iter.ci>=maxit)) {
					warning("Convergence in estimating profile likelihood CI or p-values not attained for all variables.\nConsider re-run with smaller maxstep and larger maxit.\n")
		}
		if(any(fit$iter.ci[,3]==-9)) warning("Numerical error in computing penalized likelihood ratio test for some parameters.\n")
		fit$prob[fit$iter.ci[,3]==-9]<-NA
				
  } else {
    fit$method.ci <- "Wald"
    fit$ci.lower <- exp(coefs + qnorm(alpha/2) * vars^0.5)
    fit$ci.upper <- exp(coefs + qnorm(1 - alpha/2) * vars^0.5)
    fit$prob <- 1 - pchisq((coefs^2/vars), 1)
  }
  names(fit$prob) <- names(fit$ci.upper) <- names(fit$ci.lower) <- cov.name
  attr(fit, "class") <- c("coxphf", "coxph")
  fit
}
