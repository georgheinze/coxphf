\name{coxphf}
\alias{coxphf}
\alias{Surv}
\title{Cox Regression with Firth's Penalized Likelihood}
\description{
  Implements Firth's penalized maximum likelihood bias reduction method  for Cox regression
  which has been shown to provide a solution in case of monotone likelihood (nonconvergence of likelihood function).
  The program fits profile penalized likelihood confidence intervals which were proved to outperform
  Wald confidence intervals.
}
\usage{
coxphf(formula = attr(data, "formula"), data = sys.parent(), 
       pl = TRUE, alpha = 0.05, maxit = 50, maxhs = 5, 
       epsilon = 1e-06, gconv=0.0001, maxstep = 0.5, firth = TRUE, adapt=NULL, 
	   penalty=0.5)       
}
\arguments{
  \item{formula}{a formula object, with the response on the left of the  operator, and the
    model terms on the right. The response must be a survival object as returned by the 'Surv' function (see its documentation in the survival package).
    }
  \item{data}{a data.frame in which to interpret the variables named in the 'formula' argument. }
  \item{pl}{specifies if confidence intervals and tests should be based on the profile penalized     log likelihood (\code{pl=TRUE}, the default) or on the Wald method (\code{pl=FALSE}).}
  \item{alpha}{the significance level (1-\eqn{\alpha} = the confidence level), 0.05 as default.}
  \item{maxit}{maximum number of iterations (default value is 50)}
  \item{maxhs}{maximum number of step-halvings per iterations (default value is 5). 
     The increments of the parameter vector in one Newton-Rhaphson iteration step are halved, 
     unless the new likelihood is greater than the old one, maximally doing \code{maxhs} halvings.}
  \item{epsilon}{specifies the maximum allowed change in standardized parameter estimates to
    declare convergence. Default value is 1e-6.}
  \item{gconv}{specifies the maximum allowed absolute value of first derivative of likelihood to
    declare convergence. Default value is 0.0001.}
  \item{maxstep}{specifies the maximum change of (standardized) parameter values allowed
    in one iteration. Default value is 0.5.}
  \item{firth}{use of Firth's penalized maximum likelihood (\code{firth=TRUE}, default) or the
    standard maximum likelihood method (\code{firth=FALSE}) for fitting the Cox model.}
  \item{adapt}{optional: specifies a vector of 1s and 0s, where 0 means that the corresponding parameter is fixed at 0, while 1 enables
     parameter estimation for that parameter. The length of adapt must be equal to the number of parameters to be estimated.}
  \item{penalty}{strength of Firth-type penalty. Defaults to 0.5.}

}
\details{
    The phenomenon of monotone likelihood in a sample causes parameter estimates of a Cox model to diverge,
    with infinite standard errors. Therefore, classical maximum likelihood analysis fails; the usual Wald confidence 
    intervals cover the whole range of real numbers. Monotone likelihood appears if there is single covariate
    or a linear combination of covariates such that at each event time, out of all individuals being at risk at that time,
    the individual with the highest (or at each event
    time the individual with the lowest) value for that covariate or linear combination experiences the event. It was shown that
    analysis by Firth's penalized likelihood method, particularly in conjunction with the computation
    of profile likelihood confidence intervals and penalized likelihood ratio tests is superior to maximum
    likelihood analysis. It completely removes the convergence problem mentioned in the paragraph on CONVERGENCE of
    the description of the function \code{coxph}. The \code{formula} may involve time-dependent effects or 
    time-dependent covariates. The response may
    be given in counting process style, but it cannot be used for multivariate failure times, as the program has no option
    to fit a robust covariance matrix. The user is responsible for the independency of observations within each risk set, i.e., 
    the same individual should not appear twice within the same risk set.
    
    The package coxphf provides a comprehensive tool to facilitate the
application of Firth's penalized likelihood method to Cox regression analysis. The core routines are written in Fortran 90,
(and to our knowledge this is the first package written in Fortran 90). Some description of the problem of monotone likelihood
and Firth's penalized likelihood method as a solution can be found the web page
\url{http://cemsiis.meduniwien.ac.at/en/kb/science-research/software/statistical-software/fccoxphf/}.
    
    Version 1.13  now includes a convergence check and issues a warning in case of non-convergence. Profile likelihood confidence intervals or 
	the estimation of the penalized likelihood ratio $p$-values can be vulnerable
	to non-convergence for numerical issues. In case of non-convergence problems, we suggest to first compare the output values iter.ci with the input parameter maxit. 
	Then, set maxstep to a smaller value, e.g., 0.1 and increase the number of allowed iterations to e.g. 500. This setting may slow down convergence for some 
	of the confidence limits, but proved robust also in extreme data sets. 
	
	}
\value{
 \item{coefficients}{the parameter estimates}
 \item{alpha}{the significance level = 1 - confidence level}
 \item{var}{the estimated covariance matrix}
 \item{df}{the degrees of freedom}
 \item{loglik}{the null and maximimized (penalized) log likelihood}
 \item{method.ties}{the ties handling method}
 \item{iter}{the number of iterations needed to converge}
 \item{n}{the number of observations}
 \item{y}{the response}
 \item{formula}{the model formula}
 \item{means}{the means of the covariates}
 \item{linear.predictors}{the linear predictors}
 \item{method}{the estimation method (Standard ML or Penalized ML)}
 \item{method.ci}{the confidence interval estimation method (Profile Likelihood or Wald)}
 \item{ci.lower}{the lower confidence limits}
 \item{ci.upper}{the upper confidence limits}
 \item{prob}{the p-values}
 \item{call}{the function call}
 \item{iter.ci}{the numbers of iterations needed for profile likelihood confidence interval estimation, and for maximizing the 
  restricted likelihood for p-value computation.}
}
\references{
Firth D (1993). Bias reduction of maximum likelihood estimates. \emph{Biometrika} 
  80:27--38.

Heinze G and Schemper M (2001). A Solution to the Problem of Monotone Likelihood in Cox Regression. \emph{Biometrics}
 57(1):114--119. 

Heinze G (1999). Technical Report 10/1999: The application of Firth's procedure to Cox and logistic regression. Section of Clinical Biometrics, Department of Medical Computer Sciences, University of Vienna, Vienna.
\url{http://cemsiis.meduniwien.ac.at/en/kb/science-research/software/statistical-software/fccoxphf/}
 
Heinze G and Ploner M (2002). SAS and SPLUS programs to perform Cox regression without convergence problems. \emph{Computer Methods
and Programs in Biomedicine 67:217-223} 

}
\author{Georg Heinze and Meinhard Ploner}
\note{
There exists an earlier version of coxphf for S-Plus, which is not able to involve time-dependent effects or the counting-process
representation of survival times.
}
\seealso{coxphfplot, coxphftest}
\examples{
# fixed covariate and monotone likelihood
time<-c(1,2,3)
cens<-c(1,1,1)
x<-c(1,1,0)
sim<-cbind(time,cens,x)
sim<-data.frame(sim)
coxphf(sim, formula=Surv(time,cens)~x) #convergence attained!
#coxph(sim, formula=Surv(time,cens)~x)  #no convergence!

# time-dependent covariate
test2 <- data.frame(list(start=c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
                stop =c(2, 3, 6, 7, 8, 9, 9, 9,14,17),
                event=c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
                x    =c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0) ))

summary( coxphf( formula=Surv(start, stop, event) ~ x, pl=FALSE, data=test2))


# time-dependent effect
# the coxphf function can handle interactions of a (fixed or time-dependent)
# covariate with time
# such that the hazard ratio can be expressed as a function of time

summary(coxphf(formula=Surv(start, stop, event)~x+x:log(stop), data=test2, pl=FALSE, firth=TRUE))

# note that coxph would treat x:log(stop) as a fixed covariate
# (computed before the iteration process)
# coxphf treats x:log(stop) as a time-dependent covariate which changes (
# for the same individual!) over time


# time-dependent effect with monotone likelihood

test3 <- data.frame(list(start=c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
                         stop =c(2, 3, 6, 7, 8, 9, 9, 9,14,17),
                         event=c(1, 0, 0, 1, 0, 1, 1, 0, 0, 0),
                         x    =c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0) ))
                         
summary( coxphf( formula=Surv(start, stop, event) ~ x+x:log(stop), pl=FALSE, maxit=400, data=test3))


# no convergence if option "firth" is turned off: 
# summary( coxphf(formula=Surv(start, stop, event) ~ x+x:log(stop), pl=F,
#                 data=test3, firth=FALSE)


data(breast)
fit.breast<-coxphf(data=breast, Surv(TIME,CENS)~T+N+G+CD)
summary(fit.breast)
}
\keyword{survival}