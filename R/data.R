#' Breast Cancer Data Set
#' 
#' Provides the breast cancer data set as used by Heinze & Schemper, 2001. 
#' The data sets contains information on 100 breast cancer patients, including:
#' survival time, survival status, Tumor stage, Nodal status, Grading, Cathepsin-D tumorexpression
#' 
#' @format  A data frame with 100 observations on the following 6 variables.
#' \describe{
#'   \item{\code{T}}{a numeric vector}
#'   \item{\code{N}}{a numeric vector}
#'   \item{\code{G}}{a numeric vector}
#'   \item{\code{CD}}{a numeric vector}
#'   \item{\code{TIME}}{a numeric vector}
#'   \item{\code{CENS}}{a numeric vector}
#'   \item{\code{pat}}{a numeric vector}
#' }
#' @references Heinze, G., and Schemper, M. 2001. A solution to the problem of monotone likelihood. Biometrics 57(1) pp. 114-119.
#' @keywords datasets
"breast"