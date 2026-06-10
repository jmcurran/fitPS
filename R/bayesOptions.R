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

logitPi = function(pi) {
  pi = unname(pi)

  if (!is.numeric(pi) || any(!is.finite(pi)) || any(pi <= 0 | pi >= 1)) {
    stop("pi must be numeric and strictly between 0 and 1")
  }

  unname(log(pi / (1 - pi)))
}

invLogitPi = function(eta) {
  eta = unname(eta)

  if (!is.numeric(eta) || any(!is.finite(eta))) {
    stop("eta must be finite and numeric")
  }

  unname(1 / (1 + exp(-eta)))
}

shapeToTau = function(shape) {
  shape = unname(shape)

  validateZetaShape(shape, "shape")
  unname(log(shape - 1))
}

tauToShape = function(tau) {
  tau = unname(tau)

  if (!is.numeric(tau) || any(!is.finite(tau))) {
    stop("tau must be finite and numeric")
  }

  unname(1 + exp(tau))
}

zizThetaToWorking = function(theta) {
  theta = unname(theta)

  if (!is.numeric(theta) || length(theta) != 2L) {
    stop("theta must be a numeric vector of length two")
  }

  c(eta = logitPi(theta[1]), tau = shapeToTau(theta[2]))
}

zizWorkingToTheta = function(working) {
  working = unname(working)

  if (!is.numeric(working) || length(working) != 2L) {
    stop("working must be a numeric vector of length two")
  }

  c(pi = invLogitPi(working[1]), shape = tauToShape(working[2]))
}

zizWorkingLogJacobian = function(working) {
  working = unname(working)
  theta = zizWorkingToTheta(working)

  unname(log(theta["pi"]) + log1p(-theta["pi"]) + working[2])
}

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
