#' Normalise Bayesian fitting options and resolve the posterior method and prior.
#'
#' @param bayesOptions A list of Bayesian fitting options, or `NULL`.
#' @param prior A prior object created by `makePrior()`.
#' @param allowedPosteriorMethods Character vector of posterior methods allowed in the current context.
#' @param defaultPosteriorMethod Default posterior method used when none is supplied.
#' @return A list containing the normalised `posteriorMethod` and validated `prior`.
#' @keywords internal
#' @noRd
normaliseBayesOptions = function(bayesOptions = NULL,
                                 prior,
                                 allowedPosteriorMethods = c("numerical", "mcmc", "laplace", "importance"),
                                 defaultPosteriorMethod = "numerical") {
  if (is.null(bayesOptions)) {
    bayesOptions = list()
  }

  if (!is.list(bayesOptions)) {
    stop("bayesOptions must be a list")
  }

  if ("method" %in% names(bayesOptions)) {
    stop("Use bayesOptions$posteriorMethod rather than bayesOptions$method")
  }

  unexpectedNames = setdiff(names(bayesOptions), c("posteriorMethod", "prior"))
  if (length(unexpectedNames) > 0L) {
    stop(
      "Unsupported bayesOptions element(s): ",
      paste(unexpectedNames, collapse = ", ")
    )
  }

  posteriorMethod = bayesOptions$posteriorMethod
  if (is.null(posteriorMethod)) {
    posteriorMethod = defaultPosteriorMethod
  }

  posteriorMethod = match.arg(posteriorMethod, allowedPosteriorMethods)

  if (!is.null(bayesOptions$prior) && !missing(prior)) {
    stop("Specify the Bayesian prior either as prior or bayesOptions$prior, not both")
  }

  priorObject = if (!is.null(bayesOptions$prior)) {
    bayesOptions$prior
  } else if (!missing(prior)) {
    prior
  } else {
    makePrior()
  }

  validateBayesPrior(priorObject)

  list(
    posteriorMethod = posteriorMethod,
    prior = priorObject
  )
}

#' Validate the internal structure of a Bayesian prior object.
#'
#' @param prior A prior object created by `makePrior()`.
#' @return `TRUE` invisibly when validation succeeds; otherwise an error is raised.
#' @keywords internal
#' @noRd
validateBayesPrior = function(prior) {
  if (!is.list(prior)) {
    stop("Bayesian prior must be a list-like psPrior object")
  }

  if (is.null(prior$range) || is.null(prior$logd)) {
    stop("Bayesian prior must contain range and logd elements")
  }

  validatePriorRange(prior$range)

  if (!is.function(prior$logd)) {
    stop("Bayesian prior logd element must be a function")
  }

  invisible(TRUE)
}

#' Transform a zero-inflation probability to the unconstrained logit scale.
#'
#' @param pi Zero-inflation probability on the natural scale.
#' @return Numeric values on the logit scale.
#' @keywords internal
#' @noRd
logitPi = function(pi) {
  pi = unname(pi)

  if (!is.numeric(pi) || any(!is.finite(pi)) || any(pi <= 0 | pi >= 1)) {
    stop("pi must be numeric and strictly between 0 and 1")
  }

  # The logit removes the bounded (0, 1) constraint so optimisation and
  # Gaussian approximations can operate on an unconstrained coordinate.
  unname(log(pi / (1 - pi)))
}

#' Transform a logit-scale zero-inflation parameter back to probability scale.
#'
#' @param eta Zero-inflation parameter on the logit working scale.
#' @return Numeric probabilities strictly between zero and one.
#' @keywords internal
#' @noRd
invLogitPi = function(eta) {
  eta = unname(eta)

  if (!is.numeric(eta) || any(!is.finite(eta))) {
    stop("eta must be finite and numeric")
  }

  unname(1 / (1 + exp(-eta)))
}

#' Transform a zeta shape parameter greater than one to an unconstrained working scale.
#'
#' @param shape Zeta shape parameter on the fitPS scale.
#' @return Numeric working-scale values.
#' @keywords internal
#' @noRd
shapeToTau = function(shape) {
  shape = unname(shape)

  validateZetaShape(shape, "shape")
  # Subtracting one before logging maps the required shape > 1 support to
  # the full real line while keeping the inverse transformation simple.
  unname(log(shape - 1))
}

#' Transform an unconstrained working parameter to the zeta shape scale.
#'
#' @param tau Zeta shape parameter on the unconstrained log(shape - 1) working scale.
#' @return Numeric zeta shape values greater than one.
#' @keywords internal
#' @noRd
tauToShape = function(tau) {
  tau = unname(tau)

  if (!is.numeric(tau) || any(!is.finite(tau))) {
    stop("tau must be finite and numeric")
  }

  unname(1 + exp(tau))
}

#' Transform natural zero-inflated zeta parameters to unconstrained working coordinates.
#'
#' @param theta Numeric vector containing natural ZIZ parameters `(pi, shape)`.
#' @return Named numeric vector `(eta, tau)`.
#' @keywords internal
#' @noRd
zizThetaToWorking = function(theta) {
  theta = unname(theta)

  if (!is.numeric(theta) || length(theta) != 2L) {
    stop("theta must be a numeric vector of length two")
  }

  c(eta = logitPi(theta[1]), tau = shapeToTau(theta[2]))
}

#' Transform zero-inflated zeta working coordinates to natural parameters.
#'
#' @param working Numeric vector or matrix of working parameters `(eta, tau)`.
#' @return Named numeric vector `(pi, shape)`.
#' @keywords internal
#' @noRd
zizWorkingToTheta = function(working) {
  working = unname(working)

  if (!is.numeric(working) || length(working) != 2L) {
    stop("working must be a numeric vector of length two")
  }

  c(pi = invLogitPi(working[1]), shape = tauToShape(working[2]))
}

#' Compute the log absolute Jacobian for the zero-inflated zeta working transformation.
#'
#' @param working Numeric vector or matrix of working parameters `(eta, tau)`.
#' @return Numeric log absolute Jacobian determinant.
#' @keywords internal
#' @noRd
zizWorkingLogJacobian = function(working) {
  working = unname(working)
  theta = zizWorkingToTheta(working)

  # For pi = logistic(eta), d pi / d eta = pi (1 - pi); for
  # shape = 1 + exp(tau), d shape / d tau = exp(tau). The log-Jacobian
  # is therefore log(pi) + log(1 - pi) + tau.
  unname(log(theta["pi"]) + log1p(-theta["pi"]) + working[2])
}

#' Normalise current and legacy Bayesian method selections.
#'
#' @param method Character string identifying a fitting or posterior method.
#' @param bayesOptions A list of Bayesian fitting options, or `NULL`.
#' @return A list containing normalised method and Bayesian options.
#' @keywords internal
#' @noRd
normaliseBayesMethod = function(method, bayesOptions = NULL) {
  method = match.arg(
    method,
    c("mle", "bayes", "integrate", "numerical", "mcmc", "laplace", "importance")
  )

  if (method %in% c("mle", "bayes")) {
    return(list(method = method, bayesOptions = bayesOptions))
  }

  posteriorMethod = if (method == "integrate") {
    "numerical"
  } else {
    method
  }

  if (is.null(bayesOptions)) {
    bayesOptions = list()
  }

  if (!is.list(bayesOptions)) {
    stop("bayesOptions must be a list")
  }

  if (!is.null(bayesOptions$posteriorMethod) && bayesOptions$posteriorMethod != posteriorMethod) {
    stop(
      "Legacy method = ",
      sQuote(method),
      " conflicts with bayesOptions$posteriorMethod = ",
      sQuote(bayesOptions$posteriorMethod)
    )
  }

  bayesOptions$posteriorMethod = posteriorMethod

  warning(
    "method = ",
    sQuote(method),
    ' is deprecated; use method = "bayes" with bayesOptions$posteriorMethod = ',
    sQuote(posteriorMethod),
    call. = FALSE
  )

  list(method = "bayes", bayesOptions = bayesOptions)
}

#' Normalise Bayesian options for an externally defined model.
#'
#' External models own their prior representation, so this helper deliberately
#' does not impose the legacy `psPrior` structure. Stage 9 currently exposes
#' MCMC remains the primary generic Bayesian engine. One-dimensional external
#' models may also opt into the generic numerical engine; Laplace and importance
#' remain secondary built-in capabilities.
#'
#' @param bayesOptions Optional list containing `posteriorMethod` and/or `prior`.
#' @param prior Optional model-specific prior supplied directly to [fit()].
#' @return A list containing `posteriorMethod` and `prior`.
#' @keywords internal
#' @noRd
normaliseExternalBayesOptions = function(bayesOptions = NULL, prior) {
  if (is.null(bayesOptions)) {
    bayesOptions = list()
  }
  if (!is.list(bayesOptions)) {
    stop("bayesOptions must be a list", call. = FALSE)
  }
  if ("method" %in% names(bayesOptions)) {
    stop("Use bayesOptions$posteriorMethod rather than bayesOptions$method", call. = FALSE)
  }

  unexpectedNames = setdiff(names(bayesOptions), c("posteriorMethod", "prior"))
  if (length(unexpectedNames) > 0L) {
    stop(
      "Unsupported bayesOptions element(s): ",
      paste(unexpectedNames, collapse = ", "),
      call. = FALSE
    )
  }

  posteriorMethod = bayesOptions$posteriorMethod
  if (is.null(posteriorMethod)) {
    posteriorMethod = "mcmc"
  }
  posteriorMethod = match.arg(posteriorMethod, c("mcmc", "numerical"))

  if (!is.null(bayesOptions$prior) && !missing(prior)) {
    stop("Specify the Bayesian prior either as prior or bayesOptions$prior, not both", call. = FALSE)
  }

  priorObject = if (!is.null(bayesOptions$prior)) {
    bayesOptions$prior
  } else if (!missing(prior)) {
    prior
  } else {
    stop(
      "external Bayesian models require an explicit model-specific prior",
      call. = FALSE
    )
  }

  list(
    posteriorMethod = posteriorMethod,
    prior = priorObject
  )
}
