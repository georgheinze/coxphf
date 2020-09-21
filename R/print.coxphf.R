print.coxphf <-
  function(
    x,   	# object of class coxphf
    ...		# dummy
  )
    ### MP and GH
    ### 2006-10
  {
    print(x$call)
    cat("Model fitted by", x$method)
    cat("\nConfidence intervals and p-values by", x$method.ci, "\n\n")
    
    out <- cbind(x$coefficients, diag(x$var)^0.5, exp(x$coefficients), 
                 x$ci.lower, x$ci.upper, qchisq(1 - x$prob, 1), x$prob)
    dimnames(out) <- list(names(x$coefficients), c("coef", "se(coef)", 
                                                   "exp(coef)", paste(c("lower", "upper"), 1 - x$alpha), "z", "p"))	
    
    if (x$method.ci != "Wald") 
      dimnames(out)[[2]][6] <- "Chisq"
    print(out)
    
    LL <- 2 * diff(x$loglik)
    cat("\nLikelihood ratio test=", LL, " on ", x$df, 
        " df, p=", 1 - pchisq(LL, x$df), ", n=", x$n, "\n\n", sep = "")
    
    invisible(x)
  }