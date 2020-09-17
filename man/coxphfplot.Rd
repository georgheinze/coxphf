\name{coxphfplot}
\alias{coxphfplot}
\title{Plot the Penalized Profile Likelhood Function}
\description{
  Plots the penalized profile likelihood for a specified parameter.
}
\usage{
coxphfplot(formula = attr(data, "formula"), data = sys.parent(), 
           profile, pitch = 0.05, limits, alpha = 0.05, 
           maxit = 50, maxhs = 5, epsilon = 1e-06, maxstep = 0.5, 
           firth = TRUE,  penalty=0.5, adapt=NULL, legend = "center", ...)
}
\arguments{
  \item{formula}{a formula object, with the response on the left of the  operator, and the
    model terms on the right. The response must be a survival object as returned by the 'Surv' function.}
  \item{data}{a data.frame in which to interpret the variables named in the 'formula' argument.}
  \item{profile}{a righthand formula specifying the plotted parameter, interaction or
    general term, e.g. \code{~ A} or \code{~ A : C}.}
  \item{pitch}{distances between the interpolated points in standard errors of
    the parameter estimate, the default value is 0.05.}
  \item{limits}{the range of the x-axis in terms of standard errors from the parameter estimate. The default values
    are the extremes of both confidence intervals, Wald and PL, plus or minus half a
    standard error, respectively.}
  \item{alpha}{the significance level (1-\eqn{\alpha} the confidence level,
    0.05 as default).}
  \item{maxit}{maximum number of iterations (default value is 50)}
  \item{maxhs}{maximum number of step-halvings per iterations (default value is 5). 
     The increments of the parameter vector in one Newton-Rhaphson iteration step are halved, 
     unless the new likelihood is greater than the old one, maximally doing \code{maxhs} halvings.}
  \item{epsilon}{specifies the maximum allowed change in penalized log likelihood to
    declare convergence. Default value is 0.0001.}
  \item{maxstep}{specifies the maximum change of (standardized) parameter values allowed
    in one iteration. Default value is 2.5.}
  \item{firth}{use of Firth's penalized maximum likelihood (\code{firth=TRUE}, default) or the
    standard maximum likelihood method (\code{firth=FALSE}) for fitting the Cox model.}
  \item{adapt}{optional: specifies a vector of 1s and 0s, where 0 means that the corresponding parameter is fixed at 0, while 1 enables
     parameter estimation for that parameter. The length of adapt must be equal to the number of parameters to be estimated.}
  \item{penalty}{strength of Firth-type penalty. Defaults to 0.5.}
  \item{legend}{if FALSE, legends in the plot would be omitted
    (default is TRUE).}
  \item{...}{other parameters to legend}
}
\details{
This function plots the profile (penalized) log likelihood of the specified parameter.
 A symmetric shape of the profile (penalized)
log likelihood (PPL) function allows use of Wald intervals, while an
asymmetric shape demands profile (penalized) likelihood intervals (Heinze
& Schemper (2001)).
}
\value{
 A matrix of dimension \eqn{m \times 3}, with \eqn{m = 1/\code{pitch} + 1}. With the default settings, \eqn{m=101}.
 The column headers are:
 \item{std}{the distance from the parameter estimate in standard errors}
 \item{x}{the parameter value}
 \item{log-likelihood}{the profile likelihood at \code{x}}
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
\seealso{coxphf}
\examples{
time<-c(1,2,3)
cens<-c(1,1,1)
x<-c(1,1,0)
sim<-cbind(time,cens,x)
sim<-data.frame(sim)
profplot<-coxphfplot(sim, formula=Surv(time,cens)~x, profile=~x)
}
\keyword{survival}

