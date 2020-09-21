#' Penalized Likelihood Ratio Test in Cox Regression
#' 
#' Performs a penalized likelihood ratio test for hypotheses within a Cox regression analysis using Firth's penalized likelihood.
#'
#'This function performs a penalized likelihood ratio test on some (or all) selected parameters.  
#'It can be used to test contrasts of parameters, or factors that are coded in dummy variables. 
#'The resulting object is of the class coxphftest and includes the information printed by the proper print method. 
#'
#' @param formula a formula object, with the response on the left of the  operator, and the model terms on the right. The response must be a survival object as returned by the 'Surv' function.
#' @param data a data.frame in which to interpret the variables named in the 'formula' argument.
#' @param test righthand formula of parameters to test (e.g. \code{~ B + D}). As default the null hypothesis that all parameters are 0 is tested.
#' @param values null hypothesis values, default values are 0. For testing the hypothesis H0: B1=1 and B4=2 and B5=0, specify \code{test= ~ B1 + B4 + B5} and \code{values=c(1, 2, 0)}.
#' @param maxit maximum number of iterations (default value is 50)
#' @param maxhs maximum number of step-halvings per iterations (default value is 5). The increments of the parameter vector in one Newton-Rhaphson iteration step are halved, unless the new likelihood is greater than the old one, maximally doing \code{maxhs} halvings.
#' @param epsilon specifies the maximum allowed change in penalized log likelihood todeclare convergence. Default value is 0.0001.
#' @param maxstep specifies the maximum change of (standardized) parameter values allowed in one iteration. Default value is 2.5.
#' @param firth use of Firth's penalized maximum likelihood (\code{firth=TRUE}, default) or the standard maximum likelihood method (\code{firth=FALSE}) for fitting the Cox model.
#' @param adapt optional: specifies a vector of 1s and 0s, where 0 means that the corresponding parameter is fixed at 0, while 1 enables 
#' parameter estimation for that parameter. The length of adapt must be equal to the number of parameters to be estimated.
#' @param penalty strength of Firth-type penalty. Defaults to 0.5.
#'
#' @return
#' \item{testcov}{the names of the tested model terms}
#' \item{loglik}{the restricted and unrestricted maximized (penalized) log likelihood}
#' \item{df}{the number of degrees of freedom related to the test}
#' \item{prob}{the p-value}
#' \item{call}{the function call}
#' \item{method}{the estimation method (penalized ML or ML)}
#' 
#' @export
#'
#' @examples
#' library(survival)
#' testdata <- data.frame(list(start=c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
#' stop =c(2, 3, 6, 7, 8, 9, 9, 9,14,17),
#' event=c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
#' x1    =c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0),
#' x2    =c(0, 1, 1, 1, 0, 0, 1, 0, 1, 0),
#' x3    =c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0)))
#' 
#' summary( coxphf( formula=Surv(start, stop, event) ~ x1+x2+x3, data=testdata))
#' 
#' # testing H0: x1=0, x2=0
#' 
#' coxphftest( formula=Surv(start, stop, event) ~ x1+x2+x3, test=~x1+x2,  data=testdata)
#' 
#' @references Firth D (1993). Bias reduction of maximum likelihood estimates. \emph{Biometrika} 80:27--38.
#' 
#' Heinze G and Schemper M (2001). A Solution to the Problem of Monotone Likelihood in Cox Regression. \emph{Biometrics} 57(1):114--119. 
#' 
#' Heinze G (1999). Technical Report 10/1999: The application of Firth's procedure to Cox and logistic regression. Section of Clinical Biometrics, Department of Medical Computer Sciences, University of Vienna, Vienna.
coxphftest <-
  function(
    formula, 
    data, 
    test=~.,    # righthand formula fuer zu testende Faktoren (z.B. ~ B+D )
    values,     # fixierte Werte fuer die Faktoren, default=0 (z.B. c(2,0))
    maxit=50, 
    maxhs=5, 
    epsilon=1e-6, 
    maxstep=0.5, 
    firth=TRUE,
    adapt=NULL,
    penalty=0.5
  ) {
    ### MP, GH, 2006-10
    ###
    
    call <- match.call()
    mf <- match.call(expand.dots =FALSE)
    m <- match(c("formula","data"), names(mf), 0L)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    
    n <- nrow(data)
    
    obj <- decomposeSurv(formula, data, sort = TRUE)
    NTDE <- obj$NTDE
    mmm <- cbind(obj$mm1, obj$timedata)
    cov.name <- obj$covnames
    k <- ncol(obj$mm1)
    ones <- matrix(1, n, k + NTDE)
    
    ## standardisierung
    sd1 <- apply(as.matrix(obj$mm1),2,sd)
    sd2 <- apply(as.matrix(obj$timedata),2,sd)
    Z.sd <- c(sd1, sd2 * sd1[obj$timeind])
    obj$mm1 <- scale(obj$mm1, FALSE, sd1)
    obj$timedata <- scale(obj$timedata, FALSE, sd2)
    mmm <- cbind(obj$mm1, obj$timedata)
    
    CARDS <- cbind(obj$mm1, obj$resp, ones, obj$timedata)   
    PARMS <- c(n, k, firth, maxit, maxhs, maxstep, epsilon, 1, 0.0001, 0, 0, 0, 0, NTDE, penalty)
    IOARRAY <- rbind(rep(1, k+NTDE), matrix(0, 2+k+NTDE, k + NTDE))
    if(!is.null(adapt)) IOARRAY[1,]<-adapt      
    if(NTDE>0)
      IOARRAY[4,(k+1):(k+NTDE)] <- obj$timeind
    storage.mode(CARDS) <- "double"
    storage.mode(PARMS) <- "double"
    storage.mode(IOARRAY) <- "double" #
    # --------------- Aufruf Fortran - Makro FIRTHCOX ----------------------------------
    value <- .Fortran("firthcox",
                      CARDS,
                      outpar = PARMS,
                      outtab = IOARRAY, PACKAGE="coxphf")
    if(value$outpar[8])
      warning("Error in routine FIRTHCOX; parms8 <> 0")
    loglik <- c(NA, value$outpar[11])
    
    ############## now run test formula ################
    formula2 <- as.formula(paste(as.character(formula)[2], 
                                 as.character(test)[2], sep="~"))
    obj2 <- decomposeSurv(formula2, data, sort = TRUE)
    
    cov.name2 <- obj2$covnames[obj2$ind]	 # Labels der Test-Fakt.
    k2 <- length(cov.name2)		# Anzahl Faktoren
    if(!all(cov.name2 %in% cov.name))
      stop("formula test is not a subset of whole formula!")
    pos <- match(cov.name2, cov.name) ## Position der Testfakt.
    IOARRAY[1, pos] <- 0
    if(!missing(values))
      IOARRAY[2, pos] <- values * Z.sd[pos]
    # --------------- Aufruf Fortran - Makro FIRTHCOX ----------------------------------
    value <- .Fortran("firthcox",
                      CARDS,
                      outpar = PARMS,
                      outtab = IOARRAY, PACKAGE="coxphf")
    if(value$outpar[8])
      warning("Error in routine FIRTHCOX; parms8 <> 0")
    loglik[1] <- value$outpar[11]
    
    testcov <- IOARRAY[2, ]
    testcov[-pos] <- NA
    names(testcov) <- cov.name
    
    ## return object of class "coxphftest"
    fit <- list(testcov = testcov, loglik = loglik, df = k2, 
                prob = 1 - pchisq(2 * diff(loglik), k2),call = match.call())
    if(firth)
      fit$method <- "Penalized ML"
    else fit$method <- "Standard ML"
    attr(fit, "class") <- "coxphftest"
    fit
  }