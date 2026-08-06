#' S3 fitted method for an object of class \code{psFit}
#'
#' @param object an object of class \code{psFit}, usually from \code{\link{fitDist}}
#' or \code{\link{fitZIDist}}.
#' @param n This parameter is \code{NULL} by default. If it is not \code{NULL}
#' then it must be either the number of fitted terms to be returned, or a vector
#' containing the desired fitted values.
#' @param type The fitted-value definition. `"plugIn"` preserves the existing
#' behaviour and evaluates probabilities at fitted parameter values.
#' `"posteriorMean"` returns posterior mean probabilities for Bayesian fits.
#' @param ... other arguments passed to \code{fitted}---not used.
#'
#' @return a named vector of fitted probabilities
#' @importFrom stats fitted
#' @export
fitted.psFit = function(object,
                         n = NULL,
                         type = c("plugIn", "posteriorMean"),
                         ...) {
  type = match.arg(type)

  if (identical(type, "posteriorMean")) {
    if (!identical(object$method, "bayes") ||
        !inherits(object$posterior, "psPosterior")) {
      stop("posteriorMean fitted values require a Bayesian psFit object")
    }

    probabilities = posteriorProbs(object, n = n)
    estimates = probabilities$estimate
    names(estimates) = probabilities$term
    return(estimates)
  }

  if (is.null(n)) {
    return(object$fitted)
  }

  if (any(n < 0)) {
    stop("The value(s) of n must be >= 0")
  }

  if (any(abs(n - floor(n)) > .Machine$double.eps)) {
    n = floor(n)
    warning("n must be integer valued, and so has been rounded down")
  }

  surveyType = object$psData$type

  if (surveyType == "S" && any(n == 0)) {
    stop("n must be >= 1 for S terms")
  }

  probabilityFunction = probfun(object)

  if (length(n) == 1L) {
    start = ifelse(surveyType == "P", 0, 1)
    end = n - ifelse(surveyType == "P", 1, 0)
    result = probabilityFunction(start:end)
  } else {
    result = probabilityFunction(sort(n))
  }

  result
}
