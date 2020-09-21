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