#' @exportS3Method print coxphftest
print.coxphftest <-
  function(
    x,		# object of class coxphftest
    ...		# dummy
  )
    ### MP, GH, 2006-10
  {
    print(x$call)
    cat("Model fitted by", x$method, "\n\nFactors fixed as follows:\n")
    
    ## output of factors and fixing values:
    print(x$testcov)
    LL <- 2 * diff(x$loglik) 
    out <- c(x$loglik[1], x$loglik[2], LL / 2)
    names(out) <- c("Restricted model", "Full model", "difference")
    
    ## output of likelihoods:
    cat("\nLikelihoods:\n")
    print(out)
    
    ## result summary:
    cat("\nLikelihood ratio test=", LL, " on ", x$df, 
        " df, p=", x$prob, "\n", sep="")
    
    invisible(x)
  }