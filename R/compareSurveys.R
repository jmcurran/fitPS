#' Compare two surveys on the basis of their shape parameters
#'
#' @aliases compare.surveys comp.survs
#'
#' @param x either an object of class \code{psData}---see \code{\link{readData}}
#' or an object of class \code{psFit}---see \code{\link{fitDist}}.
#' @param y either an object of class \code{psData}---see \code{\link{readData}}
#' or an object of class \code{psFit}---see \code{\link{fitDist}}.
#' @param xname an optional name for the first survey object.
#' @param yname an optional name for the second survey object.
#' @param alternative one of \code{"two.sided"}, \code{"less"}, or \code{"greater"}, depending on the type of
#' hypothesis test you wish to carry out. These may be replaced by single letter (or more) abbreviations.
#' @param null.value the true value of the difference in the shape parameters under the null hypothesis.
#' @param print if \code{TRUE} then the function will print summary output to the screen. This lets output be suppressed
#' in situations where the user wants the function to run silently.
#' @param \ldots further arguments to be passed to or from methods.
#'
#' @details
#' This function **only** works for the zeta distribution. It does not work for the zero-inflated zeta distribution. If
#' the results from fitting ZIZ models are passed to this function, then it will ignore the zero-inflated part and simply refit a zeta model.
#'
#' There is very little reason for \code{null.value} to be set to be anything other than \code{0}. However it has been included for flexibility.
#'
#' \code{alternative = "greater"} is the alternative that x has a larger shape parameter than y.
#' \code{alternative = "less"} is the alternative that x has a smaller shape parameter than y.
#'
#' @return The function returns a \code{list} of class \code{"htest"} with the following elements:
#' \describe{
#'  \item{\code{statistic}}{ -- the test statistic.}
#'  \item{\code{p.value}}{ -- the P-value associated with the estimate.}
#'  \item{\code{estimate}}{ -- the estimated difference in the shape parameters.}
#'  \item{\code{null.value}}{ -- the specified hypothesized value of the difference in shape parameters---\code{0} by default.}
#'  \item{\code{stderr}}{ -- the standard error of the difference.}
#'  \item{\code{alternative}}{ -- a character string describing the alternative hypothesis.}
#'  \item{\code{method}}{ -- a character string describing the method.}
#'  \item{\code{data.name}}{ -- a character string with the names of the two input data sets separated by " and ".}
#' }
#'
#' @examples
#' data(Psurveys)
#' lau = Psurveys$lau
#' jackson = Psurveys$jackson
#' compareSurveys(lau, jackson)
#'
#' ## Example with fitted objects - note the function just refits the models
#' fit.lau = fitDist(lau)
#' fit.jackson = fitDist(jackson)
#' compareSurveys(fit.lau, fit.jackson)
#'
#' ## Example with a bigger difference
#' compareSurveys(Psurveys$roux, lau)
#'
#' @importFrom stats pnorm
#'
#' @export
compareSurveys = function(x, ...){
  UseMethod("compareSurveys", x)
}

#' @describeIn compareSurveys Compare two surveys on the basis of their shape parameters
#' @export
compareSurveys.default = function(x,
                                  y,
                                  xname = NULL,
                                  yname = NULL,
                                  alternative = c("two.sided", "less", "greater"),
                                  null.value = 0,
                                  print = TRUE, ...){
  if(!is(x, "psData") || !is(y, "psData")){
    stop("x and y must both be objects of type psData")
  }

  fit.x = fitDist(x)
  fit.y = fitDist(y)

  shape.x = fit.x$shape
  shape.y = fit.y$shape

  v.x = fit.x$var.shape
  v.y = fit.y$var.shape

  se.diff = sqrt(v.x + v.y)
  z = (shape.x - shape.y - null.value) / se.diff

  alternative = match.arg(alternative)

  P = switch(alternative,
              two.sided = 2 * (1 - pnorm(abs(z))),
              less = pnorm(z),
              greater = pnorm(z, lower.tail = FALSE)
      )

  estimate = c(shape.x, shape.y)
  names(estimate) = c(paste("Shape of", ifelse(is.null(xname), deparse1(substitute(x)), xname)),
                      paste("Shape of", ifelse(is.null(yname), deparse1(substitute(y)), yname)))
  names(null.value) = "difference in shape parameters"
  names(z) = "z"

  rval = list(statistic = z,
              p.value = P,
              estimate = estimate,
              null.value = null.value,
              stderr = se.diff,
              alternative = alternative,
              method = paste0(ifelse(alternative == "two.sided", "Two-sided", "One-sided"), " Wald test"),
              data.name = paste(ifelse(is.null(xname), deparse1(substitute(x)), xname),
                                "and",
                                ifelse(is.null(yname), deparse1(substitute(y)), yname))
  )
  class(rval) = "htest"

  if(print){
    print(rval)
  }

  invisible(rval)
}

#' @describeIn compareSurveys Compare two surveys on the basis of their shape parameters
#' @export
compareSurveys.psData = function(x, y, ...){
  compareSurveys.default(x,
                         y,
                         xname = deparse1(substitute(x)),
                         yname = deparse1(substitute(y)),
                         ...)
}

#' @describeIn compareSurveys Compare two surveys on the basis of their shape parameters
#' @export
compareSurveys.psFit = function(x, y, ...){
  compareSurveys.default(x$psData,
                         y$psData,
                         xname = paste0(deparse1(substitute(x)), "$psData"),
                         yname = paste0(deparse1(substitute(y)), "$psData"),
                         ...)
}

#' @describeIn compareSurveys Compare two surveys on the basis of their shape parameters
compare.surveys = compareSurveys

#' @describeIn compareSurveys Compare two surveys on the basis of their shape parameters
comp.survs = compareSurveys


