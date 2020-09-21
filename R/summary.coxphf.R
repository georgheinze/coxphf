#' @exportS3Method summary coxphf
summary.coxphf <-
  function(
    object,		# object of class coxphf
    ...			# dummy
  )
    ### MP and GH
    ### 2006-10
  {
    print(object$call)
    cat("\nModel fitted by", object$method)
    cat("\nConfidence intervals and p-values by", object$method.ci, "\n\n")
    
    out <- cbind(object$coefficients, diag(object$var)^0.5, 
                 exp(object$coefficients), 
                 object$ci.lower, object$ci.upper,  
                 qchisq(1 - object$prob, 1), object$prob)
    dimnames(out) <- list(names(object$coefficients), c("coef", "se(coef)", 
                                                        "exp(coef)", paste(c("lower", "upper"), 1 - object$alpha), "z", "p"))
    if (object$method.ci != "Wald") 
      dimnames(out)[[2]][6] <- "Chisq"
    
    print(out)
    
    LL <- 2 * diff(object$loglik)
    cat("\nLikelihood ratio test=", LL, " on ", object$df, 
        " df, p=", 1 - pchisq(LL, object$df), ", n=", object$n, sep = "")
    wald.z <- t(coef(object)) %*% solve(object$var) %*% coef(object) 
    cat("\nWald test =", wald.z, "on", object$df, "df, p =", 
        1 - pchisq(wald.z, object$df))
    cat("\n\nCovariance-Matrix:\n")
    print(object$var)
    
    invisible(object)
  }