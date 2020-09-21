#' @name coxphf
#' @title Cox Regression with Firth's Penalized Likelihood
#' @docType package
#' @aliases {coxphf}-package
#' 
#' @description 
#' Implements Firth's penalized maximum likelihood bias reduction method  for Cox regression which has been 
#' shown to provide a solution in case of monotone likelihood (nonconvergence of likelihood function).
#' The program fits profile penalized likelihood confidence intervals which were proved to outperform
#' Wald confidence intervals.
#' 
#' @details 
#' The phenomenon of monotone likelihood in a sample causes parameter estimates of a Cox model to diverge, 
#' with infinite standard errors. Therefore, classical maximum likelihood analysis fails; the usual Wald confidence 
#' intervals cover the whole range of real numbers. Monotone likelihood appears if there is single covariate 
#' or a linear combination of covariates such that at each event time, out of all individuals being at risk at that time, 
#' the individual with the highest (or at each event time the individual with the lowest) value for that covariate 
#' or linear combination experiences the event. It was shown that analysis by Firth's penalized likelihood method, 
#' particularly in conjunction with the computation of profile likelihood confidence intervals and penalized 
#' likelihood ratio tests is superior to maximum likelihood analysis. It completely removes the convergence
#'  problem mentioned in the paragraph on CONVERGENCE of the description of the function \code{coxph}. 
#'  The \code{formula} may involve time-dependent effects or time-dependent covariates. The response may be given 
#'  in counting process style, but it cannot be used for multivariate failure times, as the program has no option 
#'  to fit a robust covariance matrix. The user is responsible for the independency of observations within each risk set, i.e., 
#'  the same individual should not appear twice within the same risk set.
#'  
#'  The package coxphf provides a comprehensive tool to facilitate the application of Firth's penalized 
#'  likelihood method to Cox regression analysis. The core routines are written in Fortran 90, (and to our knowledge this is the first package written in Fortran 90). Some description of the problem of monotone likelihood
#'  and Firth's penalized likelihood method as a solution can be found the web page 
#'  \url{https://cemsiis.meduniwien.ac.at/en/kb/science-research/software/statistical-software/fccoxphf/}.
#'  
#'  Version 1.13  now includes a convergence check and issues a warning in case of non-convergence. Profile likelihood confidence intervals or 
#'  the estimation of the penalized likelihood ratio $p$-values can be vulnerable
#'  non-convergence for numerical issues. In case of non-convergence problems, we suggest to first compare the output values iter.ci with the input parameter maxit. 
#'  Then, set maxstep to a smaller value, e.g., 0.1 and increase the number of allowed iterations to e.g. 500. This setting may slow down convergence for some 
#'  of the confidence limits, but proved robust also in extreme data sets. 
#' 
#' @author Georg Heinze <georg.heinze@meduniwien.ac.at> and Meinhard Ploner
#' @references 
#' 
#' 
#' @keywords survival
#' 
#' @importFrom survival Surv
#' @importFrom graphics axis mtext par plot points segments title
#' @importFrom stats as.formula coef model.extract model.frame model.matrix pchisq qchisq qnorm sd model.response
#' @importFrom utils tail
#' 
#' @useDynLib coxphf, .registration=TRUE
"_PACKAGE"