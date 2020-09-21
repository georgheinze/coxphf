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