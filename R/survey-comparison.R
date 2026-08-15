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

  fit.x = fit(x, model = zetaModel())
  fit.y = fit(y, model = zetaModel())

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
#' Compare two or more surveys on the basis of their shape parameters using a Likelihood Ratio Test
#'
#' @param \ldots two or more objects of class \code{"psData"}---see \code{\link{readData}}.
#'
#' @details
#' This function **only** works for the zeta distribution. The function carries out a likelihood ratio test (LRT) to test the null hypothesis
#' \deqn{H_0: \alpha_1 = \alpha_2 = \ldots = \alpha_K}{H_0: alpha_1 = alpha_2 = \ldots = alpha_K} versus the alternative
#' \deqn{H_1: \alpha_i \neq \alpha_j \mbox{ for some } i \neq j \in \left\{1, \ldots, K\right\},}{H_1: alpha_i != alpha_j for some  i != j in {1, \ldots, K},}
#' where \eqn{\alpha_i}{alpha_i} is the shape parameter for the zeta distribution of the \eqn{i^\mathrm{th}}{ith} survey.
#'
#' @return The function returns a \code{list} of class \code{"htest"} with the following elements:
#' \describe{
#'  \item{\code{statistic}}{ -- the test statistic.}
#'  \item{\code{parameter}}{ -- the degrees of freedom for the test}
#'  \item{\code{p.value}}{ -- the P-value associated with the estimate.}
#'  \item{\code{method}}{ -- a character string describing the method hypothesis.}
#'  \item{\code{data.name}}{ -- the names of the data sets used in the test}
#' }
#'
#' @examples
#' data(Psurveys)
#' lau = Psurveys$lau
#' jackson = Psurveys$jackson
#' compareSurveysLRT(lau, jackson)
#'
#' ## Example with three surveys
#' roux = Psurveys$roux
#' compareSurveysLRT(lau, jackson, roux)
#'
#' @importFrom stats pchisq
#'
#' @export
compareSurveysLRT = function(...){
  Surveys = list(...)

  if(length(Surveys) < 2){
    stop("You need to supply at least two surveys to carry out this test")
  }

  if(any(sapply(Surveys, function(x){
    !is(x, "psData")
  }))){
    stop("All the arguments supplied to this function must be of class psData")
  }

  types = sapply(Surveys, function(x){
    x$type
  })

  if(length(unique(types)) > 1){
    stop("All the surveys must be of the same type - either P or S")
  }

  fits = lapply(Surveys, function(survey) {
    fit(survey, model = zetaModel())
  })

  ll.mle = sum(sapply(fits, function(x){
    -x$fit$value
  }))

  ll.H0 = -fit(combineSurveys(...), model = zetaModel())$fit$value

  LRT.stat = 2 * (ll.mle - ll.H0)
  names(LRT.stat) = "X-squared"
  parameter = length(Surveys) - 1
  names(parameter) = "df"

  P = pchisq(LRT.stat, df = length(Surveys) - 1, lower.tail = FALSE)
  dname = paste(unlist(match.call(expand.dots = FALSE)$...), collapse = ", ")

  rval = list(statistic = LRT.stat,
              parameter = length(Surveys) - 1,
              p.value = P,
              method = "Likelihood Ratio Test",
              data.name = dname
  )
  class(rval) = "htest"
  return(rval)
}
#' Fit zeta and zero-inflated zeta models for internal comparison.
#'
#' @param x An input object or numeric vector required by the helper.
#' @param start Starting values for optimisation or fitting.
#' @param ... Additional arguments passed to the underlying fitting or helper routine.
#' @return An internal comparison object containing fitted candidate models.
#' @keywords internal
#' @noRd
fitCompare = function(x, start = list(zeta = 1, ziz = c(0.5, 1)), ...){
  fit.zeta = fit(x, model = zetaModel(), start = start$zeta, ...)
  fit.ziz = fit(x, model = zizModel(), start = start$ziz, ...)
  p.zeta = probfun(fit.zeta)
  p.ziz = probfun(fit.ziz)
  raw = x$data$rn / sum(x$data$rn)

  nmax = max(x$data$n) + 1

  if(x$type == "P"){
    fitted = data.frame(raw = rep(0, nmax + 1),
                        zeta = rep(0, nmax + 1),
                        ziz = rep(0, nmax + 1))
    fitted$raw[x$data$n + 1] = raw
    fitted$zeta = p.zeta(0:nmax)
    fitted$ziz = p.ziz(0:nmax)
  }else{
    fitted = data.frame(raw = rep(0, nmax),
                        zeta = rep(0, nmax),
                        ziz = rep(0, nmax))
    fitted$raw[x$data$n] = raw
    fitted$zeta = p.zeta(1:nmax)
    fitted$ziz = p.ziz(1:nmax)
  }

  print(fitted)
  invisible(list(raw = raw, fit.zeta = fit.zeta, fit.ziz = fit.ziz, fitted = fitted))
}
