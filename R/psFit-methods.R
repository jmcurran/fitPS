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
#' `"bootstrapMean"` returns bootstrap mean probabilities when a bootstrap
#' distribution has been attached to the fit.
#' @param ... other arguments passed to \code{fitted}---not used.
#'
#' @return a named vector of fitted probabilities
#' @importFrom stats fitted
#' @export
fitted.psFit = function(object,
                         n = NULL,
                         type = c("plugIn", "posteriorMean", "bootstrapMean"),
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

  if (identical(type, "bootstrapMean")) {
    if (!inherits(object$bootstrap, "psBootstrap")) {
      stop(
        "bootstrapMean fitted values require a psFit object with an attached ",
        "bootstrap; use fit(..., method = 'bootstrap')"
      )
    }

    probabilities = bootstrapProbs(object, n = n)
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
#' Extract the maximised log-likelihood from a fitPS model fit
#'
#' Returns the model log-likelihood evaluated at the maximum-likelihood estimate.
#' The returned `logLik` object includes the model parameter count and number of
#' observations, allowing [stats::AIC()] and [stats::BIC()] to use the shared
#' fitPS model-comparison contract.
#'
#' Parametric Bayesian fits are rejected because their stored representative
#' parameter values are posterior summaries rather than maximum-likelihood
#' estimates. Bootstrap fits retain the underlying maximum-likelihood point fit
#' and therefore remain valid for AIC and BIC calculations.
#'
#' @param object An object of class `psFit`.
#' @param ... Additional arguments retained for S3 compatibility.
#' @return An object of class `logLik` with `df` and `nobs` attributes.
#' @export
logLik.psFit = function(object, ...) {
  if (!is(object, "psFit")) {
    stop("object must be an object of class psFit")
  }
  if (!object$method %in% c("mle", "bootstrap")) {
    stop("AIC and BIC require a maximum-likelihood psFit object", call. = FALSE)
  }

  model = modelFromFit(object)
  parameters = fitModelParameters(object, model)
  value = modelLogLikelihood(model, parameters = parameters, data = object$psData)
  attr(value, "df") = modelParameterCount(model)
  attr(value, "nobs") = sum(object$psData$data$rn)
  class(value) = "logLik"
  value
}
#' S3 print method for an object of class \code{psFit}
#'
#' @param x An object of class \code{psFit}, usually returned by [fit()] or one
#'   of the deprecated compatibility fitters.
#' @param nterms Number of probability terms to print. If \code{NULL}, print
#'   the terms already stored in the fitted object, capped at 10 for posterior
#'   and bootstrap summaries.
#' @param ... Other arguments passed to delegated print methods.
#'
#' @return No return value, called for side effects.
#' @export
print.psFit = function(x, nterms = NULL, ...) {
  isBayes = identical(x$method, "bayes")

  if (isBayes && inherits(x$posterior, "psPosterior")) {
    print(x$posterior, nterms = nterms, ...)
    return(invisible(x))
  }

  if (x$model == "ziz") {
    cat("The estimated mixing parameter, pi, is", signif(x$pi, 4), "\n")
  }

  if (x$model %in% c("zeta", "ziz")) {
    cat("The estimated shape parameter is", round(x$shape, 4), "\n")
  }

  if (x$model %in% c("log", "logarithmic")) {
    cat("The estimated model parameter pi is", signif(x$pi, 4), "\n")
  }

  if (x$model == "zeta") {
    cat(
      "The standard error of shape parameter is",
      round(sqrt(x$var.shape), 4),
      "\n"
    )
  }

  if (!x$model %in% c("zeta", "ziz", "log", "logarithmic")) {
    model = modelFromFit(x)
    parameters = fitModelParameters(x, model)
    cat("Estimated model parameters:\n")
    print(parameters)
  }

  if (is.null(nterms)) {
    nterms = length(x$fitted)
  }
  if (!is.numeric(nterms) || length(nterms) != 1L || !is.finite(nterms) ||
      nterms < 1 || nterms != floor(nterms)) {
    stop("nterms must be one positive integer")
  }

  fittedValues = fitted(x, n = as.integer(nterms), type = "plugIn")
  cat("The first", nterms, "fitted values are:\n")
  print(fittedValues)

  if (inherits(x$bootstrap, "psBootstrap")) {
    diagnostics = x$bootstrap$diagnostics
    cat("\nBootstrap uncertainty is attached")
    if (!is.null(diagnostics$B)) {
      cat(" (", diagnostics$B, " replicates)", sep = "")
    }
    cat(".\n")
    cat("Use bootstrapProbs(), summary(), or plot(x$bootstrap) for details.\n")
  }

  invisible(x)
}
#' S3 summary method for an object of class \code{psFit}
#'
#' @param object An object of class \code{psFit}, usually returned by [fit()] or
#'   one of the deprecated compatibility fitters.
#' @param nterms Optional number of leading posterior or bootstrap probability
#'   summaries to include when an uncertainty object is attached.
#' @param ... Other arguments passed to delegated summary methods.
#'
#' @details
#' Bayesian fits delegate to the attached \code{psPosterior} object. For
#' maximum-likelihood fits, built-in zeta, ZIZ, and logarithmic presentation is
#' preserved. Other models use the public model contract to report their fitted
#' natural-scale parameters and, when available, standard errors obtained from
#' the fitted covariance matrix. If a \code{psBootstrap} object is attached,
#' its parameter and probability summaries are printed afterwards.
#'
#' @importFrom stats pchisq printCoefmat
#' @return For Bayesian fits, a \code{summary.psPosterior} object. For
#'   maximum-likelihood fits, invisibly returns the coefficient matrix.
#' @export
summary.psFit = function(object, nterms = NULL, ...) {
  if (identical(object$method, "bayes") &&
      inherits(object$posterior, "psPosterior")) {
    return(summary(object$posterior, nterms = nterms, ...))
  }

  if (object$model == "zeta") {
    cmat = matrix(c(object$shape, sqrt(object$var.shape)), nrow = 1)
    colnames(cmat) = c("Estimate", "Std.Err")
    rownames(cmat) = "shape"
  } else if (object$model == "ziz") {
    llH0 = -fit(object$psData, model = zetaModel())$fit$value
    llMle = -object$fit$value
    lrtStat = 2 * (llMle - llH0)

    cmat = rbind(
      c(object$shape, sqrt(object$var.cov[2, 2]), NA, NA),
      c(
        object$pi,
        sqrt(object$var.cov[1, 1]),
        lrtStat,
        pchisq(lrtStat, 1, lower.tail = FALSE)
      )
    )
    colnames(cmat) = c("Estimate", "Std.Err", "X^2 value", "Pr(>X^2)")
    rownames(cmat) = c("shape", "pi")
  } else if (object$model %in% c("log", "logarithmic")) {
    cmat = matrix(c(object$pi, sqrt(object$var.pi)), nrow = 1L)
    colnames(cmat) = c("Estimate", "Std.Err")
    rownames(cmat) = "pi"
  } else {
    model = modelFromFit(object)
    estimates = fitModelParameters(object, model)
    standardErrors = rep(NA_real_, length(estimates))
    names(standardErrors) = names(estimates)

    covariance = object$var.cov
    if (is.matrix(covariance) &&
        all(names(estimates) %in% rownames(covariance)) &&
        all(names(estimates) %in% colnames(covariance))) {
      variances = diag(covariance[names(estimates), names(estimates), drop = FALSE])
      valid = is.finite(variances) & variances >= 0
      standardErrors[valid] = sqrt(variances[valid])
    }

    cmat = cbind(
      Estimate = as.numeric(estimates),
      Std.Err = as.numeric(standardErrors)
    )
    rownames(cmat) = names(estimates)
  }

  printCoefmat(cmat)

  if (inherits(object$bootstrap, "psBootstrap")) {
    cat("\n")
    print(summary(object$bootstrap, nterms = nterms, ...))
  }

  invisible(cmat)
}
