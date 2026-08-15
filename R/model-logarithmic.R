# Logarithmic model implementation
#
# Distribution-specific logarithmic fitting, probability, posterior, and
# model-dispatch code is intentionally consolidated in this file.

# ---- fitlogDist.R ----
#' Fit a logarithmic distribution to forensic data
#'
#' Fits the logarithmic-series distribution to P- or S-survey data using
#' maximum likelihood or the common fitPS Bayesian posterior architecture.
#' The logarithmic model has probability mass function
#' \deqn{p(k) = -\frac{\pi^k}{k\log(1-\pi)}, \quad 0 < \pi < 1.}
#'
#' @param x An object of class `psData`.
#' @param nterms Number of fitted P/S probability terms to retain.
#' @param method Fitting method, either `"mle"` or `"bayes"`.
#' @param prior Optional `psPrior` used for Bayesian fitting. When omitted, a
#'   uniform prior on `(0.001, 0.999)` is used.
#' @param bayesOptions Optional Bayesian controls. The logarithmic model supports
#'   `"numerical"` and `"mcmc"` posterior engines.
#' @param ... Additional fitting controls. For MLE, `start` is an optional
#'   initial value in `(0, 1)`. For MCMC, `pi0`, `nIter`, `nBurnIn`, and
#'   `silent` are passed to the logarithmic MCMC model method.
#' @section Deprecated:
#' `fitlogDist()` is retained as a compatibility wrapper. New code should use
#' `fit(x, model = logarithmicModel())`.
#'
#' @return An object of class `psFit`.
#' @importFrom stats optim
#' @importFrom VGAM dlog
#' @export
#'
#' @examples
#' data(Psurveys)
#' fit = fitlogDist(Psurveys$roux)
#' fit
fitlogDist = function(x,
                       nterms = 10,
                       method = c("mle", "bayes"),
                       prior,
                       bayesOptions = NULL,
                       ...) {
  signalDeprecatedFitter(
    old = "fitlogDist",
    replacement = "fit(x, model = logarithmicModel())"
  )

  method = match.arg(method)
  args = list(
    x = x,
    model = logarithmicModel(),
    nterms = nterms,
    method = method,
    bayesOptions = bayesOptions,
    ...
  )
  if (!missing(prior)) {
    args$prior = prior
  }

  do.call(fit, args)
}

#' Internal logarithmic-series fitting implementation.
#'
#' @inheritParams fitlogDist
#' @return An object of class `psFit`.
#' @keywords internal
#' @noRd
fitlogDistImpl = function(x,
                       nterms = 10,
                       method = c("mle", "bayes"),
                       prior,
                       bayesOptions = NULL,
                       ...) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData")
  }
  modelObservationData(logarithmicModel(), x)

  method = match.arg(method)
  model = logarithmicModel()
  nvals = posteriorProbabilityIndices(x$type, nterms)

  if (identical(method, "mle")) {
    dotargs = list(...)
    start = if ("start" %in% names(dotargs)) {
      dotargs$start
    } else if ("pi" %in% names(dotargs)) {
      dotargs$pi
    } else {
      0.5
    }
    validateLogarithmicPi(start, "start")

    objective = function(pi) {
      -modelLogLikelihood(model, parameters = list(pi = pi), data = x)
    }
    fit = optim(
      par = start,
      fn = objective,
      method = "L-BFGS-B",
      lower = sqrt(.Machine$double.eps),
      upper = 1 - sqrt(.Machine$double.eps),
      hessian = TRUE
    )
    pi = unname(fit$par[1L])
    fittedMatrix = modelProbabilities(
      model,
      parameters = list(pi = pi),
      n = nvals,
      type = x$type
    )
    fitted = as.numeric(fittedMatrix[1L, ])
    names(fitted) = colnames(fittedMatrix)

    result = list(
      psData = x,
      fit = fit,
      pi = pi,
      var.pi = 1 / fit$hessian[1L, 1L],
      fitted = fitted,
      model = model$model,
      modelObject = model,
      method = "mle"
    )
    class(result) = "psFit"
    return(result)
  }

  options = if (missing(prior)) {
    if (is.null(bayesOptions)) {
      bayesOptions = list()
    }
    if (is.null(bayesOptions$prior)) {
      bayesOptions$prior = makePrior(
        family = "uniform",
        range = c(0.001, 0.999)
      )
    }
    normaliseBayesOptions(bayesOptions = bayesOptions)
  } else {
    normaliseBayesOptions(bayesOptions = bayesOptions, prior = prior)
  }

  engine = posteriorEngine(options$posteriorMethod)
  validateEngineModelPair(engine, model)
  validateLogarithmicPriorRange(options$prior$range)

  result = fitBayesianModel(
    model = model,
    posteriorMethod = options$posteriorMethod,
    x = x,
    prior = options$prior,
    nterms = nterms,
    ...
  )
  result$bayesOptions = options
  result
}

#' @rdname fitlogDist
#' @export
fitLogdist = fitlogDist

#' @rdname fitlogDist
#' @export
fitlogdist = fitlogDist

# ---- logarithmicProbabilities.R ----
#' Evaluate logarithmic P/S probabilities
#'
#' This helper evaluates the logarithmic-series distribution on the latent
#' positive-integer support used by fitPS. P probabilities use `n + 1`, while
#' S probabilities use `n` directly.
#'
#' @param pi Numeric logarithmic-distribution parameter values in `(0, 1)`.
#' @param n Requested P/S probability indices.
#' @param type Survey type, either `"P"` or `"S"`.
#' @return Numeric matrix with one row per parameter value and one column per
#'   requested probability term.
#' @importFrom VGAM dlog
#' @keywords internal
#' @noRd
logarithmicProbabilities = function(pi, n, type) {
  if (!is.numeric(pi) || length(pi) == 0L || any(!is.finite(pi))) {
    stop("pi must contain finite numeric values")
  }
  if (any(pi <= 0 | pi >= 1)) {
    stop("pi must lie strictly between 0 and 1")
  }

  type = normaliseSurveyType(type)
  n = normaliseProbabilityIndices(n, type)
  support = latentPsValues(n, type)

  values = outer(
    pi,
    support,
    Vectorize(function(parameter, value) {
      dlog(value, shape = parameter)
    })
  )
  colnames(values) = psProbabilityTermNames(n, type)
  values
}

#' Validate the logarithmic-distribution parameter
#'
#' @param pi Candidate scalar parameter.
#' @param name Name used in error messages.
#' @return `pi`, invisibly, when valid.
#' @keywords internal
#' @noRd
validateLogarithmicPi = function(pi, name = "pi") {
  if (!is.numeric(pi) || length(pi) != 1L || !is.finite(pi) ||
      pi <= 0 || pi >= 1) {
    stop(name, " must be one finite numeric value strictly between 0 and 1")
  }
  invisible(pi)
}

# ---- logarithmicPosterior.R ----
#' Validate logarithmic posterior prior support
#'
#' @param range Two-element prior support.
#' @return `range`, invisibly, when valid.
#' @keywords internal
#' @noRd
validateLogarithmicPriorRange = function(range) {
  validatePriorRange(range)
  if (range[1L] <= 0 || range[2L] >= 1) {
    stop("logarithmic prior range must lie strictly inside (0, 1)")
  }
  invisible(range)
}

#' Fit a numerical posterior for the logarithmic model
#'
#' @param model A `logarithmicModel` descriptor.
#' @param engine A numerical posterior engine.
#' @param x An object of class `psData`.
#' @param prior A `psPrior` object.
#' @param ... Additional controls reserved for future use.
#' @return A `numericalPosteriorRepresentation` object.
#' @keywords internal
#' @exportS3Method fitNumericalPosteriorModel logarithmicModel
#' @noRd
fitNumericalPosteriorModel.logarithmicModel = function(model,
                                                        engine,
                                                        x,
                                                        prior,
                                                        ...) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData")
  }
  if (!inherits(prior, "psPrior")) {
    stop("prior must be an object of class psPrior")
  }

  modelObservationData(model, x)
  validateLogarithmicPriorRange(prior$range)

  lower = prior$range[1L]
  upper = prior$range[2L]
  # The prior support helper is intentionally strict at its endpoints. Inset
  # the optimisation bounds so L-BFGS-B never evaluates a log prior of -Inf
  # merely because it lands exactly on a prior boundary. Integration can still
  # use the full interval because endpoint values have measure zero.
  boundaryInset = sqrt(.Machine$double.eps) * (upper - lower)
  optimisationLower = lower + boundaryInset
  optimisationUpper = upper - boundaryInset
  negLogPosterior = function(pi) {
    -(
      modelLogLikelihood(model, parameters = list(pi = pi), data = x) +
        prior$logd(pi)
    )
  }

  optimum = optim(
    par = mean(c(lower, upper)),
    fn = negLogPosterior,
    method = "L-BFGS-B",
    lower = optimisationLower,
    upper = optimisationUpper
  )
  logShift = optimum$value

  unnormalisedDensity = Vectorize(function(pi) {
    exp(-negLogPosterior(pi) + logShift)
  })
  normalisingIntegral = integrate(
    unnormalisedDensity,
    lower = lower,
    upper = upper
  )
  logNormalisingConstant = log(normalisingIntegral$value)
  density = Vectorize(function(pi) {
    exp(-negLogPosterior(pi) + logShift - logNormalisingConstant)
  })

  meanIntegral = integrate(
    function(pi) {
      pi * density(pi)
    },
    lower = lower,
    upper = upper
  )
  secondMomentIntegral = integrate(
    function(pi) {
      pi^2 * density(pi)
    },
    lower = lower,
    upper = upper
  )

  posteriorMean = meanIntegral$value
  posteriorVariance = secondMomentIntegral$value - posteriorMean^2

  newPsPosteriorRepresentation(
    engine = engine,
    value = list(
      density = density,
      mean = c(pi = posteriorMean),
      variance = matrix(
        posteriorVariance,
        nrow = 1L,
        dimnames = list("pi", "pi")
      )
    ),
    metadata = list(
      model = model$model,
      bounds = c(lower = lower, upper = upper),
      mode = c(pi = unname(optimum$par)),
      normalisingError = normalisingIntegral$abs.error,
      meanError = meanIntegral$abs.error,
      secondMomentError = secondMomentIntegral$abs.error
    )
  )
}

#' Fit an MCMC posterior for the logarithmic model
#'
#' @param model A `logarithmicModel` descriptor.
#' @param engine An MCMC posterior engine.
#' @param x An object of class `psData`.
#' @param prior A `psPrior` object.
#' @param pi0 Initial logarithmic parameter value.
#' @param nIter Number of retained MCMC samples.
#' @param nBurnIn Number of burn-in iterations.
#' @param silent Logical; suppress progress output when `TRUE`.
#' @param ... Additional controls reserved for future use.
#' @return An `mcmcPosteriorRepresentation` object.
#' @keywords internal
#' @exportS3Method fitMcmcPosteriorModel logarithmicModel
#' @noRd
fitMcmcPosteriorModel.logarithmicModel = function(model,
                                                   engine,
                                                   x,
                                                   prior,
                                                   pi0 = 0.5,
                                                   nIter = 1e4,
                                                   nBurnIn = 1e3,
                                                   silent = TRUE,
                                                   ...) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData")
  }
  if (!inherits(prior, "psPrior")) {
    stop("prior must be an object of class psPrior")
  }

  validateLogarithmicPi(pi0, "pi0")
  validateLogarithmicPriorRange(prior$range)
  modelObservationData(model, x)

  if (nIter < 1000) {
    warning("The number of samples from the MCMC chain really should be 1000 or higher.")
  }
  if (nIter <= 0 || nBurnIn <= 0) {
    stop("nIter and nBurnIn must be greater than zero.")
  }

  lower = prior$range[1L]
  upper = prior$range[2L]
  nTotal = nIter + nBurnIn
  proposals = runif(nTotal, lower, upper)
  logUniforms = log(runif(nTotal))
  chain = numeric(nIter)

  logPosterior = function(pi) {
    modelLogLikelihood(model, parameters = list(pi = pi), data = x) +
      prior$logd(pi)
  }

  currentLogPosterior = logPosterior(pi0)
  if (!is.finite(currentLogPosterior)) {
    stop("Log likelihood is not finite at starting value")
  }

  if (!silent) {
    progress = txtProgressBar(
      min = 1,
      max = nTotal,
      initial = 1,
      style = 3,
      label = "Burning in"
    )
  }

  i = 1L
  while (i <= nTotal) {
    proposedPi = proposals[i]
    proposedLogPosterior = logPosterior(proposedPi)

    if (proposedLogPosterior > currentLogPosterior ||
        logUniforms[i] < proposedLogPosterior - currentLogPosterior) {
      pi0 = proposedPi
      currentLogPosterior = proposedLogPosterior
    }

    if (i > nBurnIn) {
      chain[i - nBurnIn] = pi0
    }

    i = i + 1L
    if (!silent) {
      if (i <= nBurnIn) {
        setTxtProgressBar(progress, i)
      } else if (i <= nTotal) {
        setTxtProgressBar(progress, i, label = "Sampling")
      }
    }
  }

  if (!silent) {
    close(progress)
  }

  posteriorMean = mean(chain)
  posteriorVariance = var(chain)
  newPsPosteriorRepresentation(
    engine = engine,
    value = list(
      chain = chain,
      mean = c(pi = posteriorMean),
      variance = matrix(
        posteriorVariance,
        nrow = 1L,
        dimnames = list("pi", "pi")
      )
    ),
    metadata = list(
      model = model$model,
      bounds = c(lower = lower, upper = upper),
      nIter = nIter,
      nBurnIn = nBurnIn,
      acceptance = mean(diff(chain) != 0)
    )
  )
}

#' Summarise numerical logarithmic posterior probabilities
#'
#' @param model A `logarithmicModel` descriptor.
#' @param engine A numerical posterior engine.
#' @param representation A numerical posterior representation.
#' @param x An object of class `psData`.
#' @param nterms Number of P/S terms to summarise.
#' @param level Equal-tailed credible interval level.
#' @param nGrid Number of quadrature-grid points used for derived probabilities.
#' @param ... Additional controls reserved for future use.
#' @return A data frame of posterior probability summaries.
#' @keywords internal
#' @exportS3Method summariseNumericalPosteriorProbabilities logarithmicModel
#' @noRd
summariseNumericalPosteriorProbabilities.logarithmicModel = function(model,
                                                                      engine,
                                                                      representation,
                                                                      x,
                                                                      nterms,
                                                                      level = 0.95,
                                                                      nGrid = 2001L,
                                                                      ...) {
  bounds = representation$metadata$bounds
  densityFunction = representation$value$density
  if (!is.numeric(bounds) || length(bounds) != 2L ||
      any(!is.finite(bounds)) || !is.function(densityFunction)) {
    stop("Numerical logarithmic posterior representation is incomplete")
  }

  nGrid = as.integer(nGrid)
  if (!is.finite(nGrid) || nGrid < 101L) {
    stop("nGrid must be at least 101")
  }

  pi = seq(bounds[1L], bounds[2L], length.out = nGrid)
  densityValues = pmax(0, as.numeric(densityFunction(pi)))
  spacing = (bounds[2L] - bounds[1L]) / (nGrid - 1L)
  weights = densityValues * spacing
  weights[c(1L, nGrid)] = weights[c(1L, nGrid)] / 2

  probabilities = modelProbabilities(
    model = model,
    parameters = list(pi = pi),
    n = posteriorProbabilityIndices(x$type, nterms),
    type = x$type
  )

  summariseZizProbabilities(
    probabilities = probabilities,
    weights = weights,
    level = level,
    posteriorMethod = posteriorEngineName(engine)
  )
}

# ---- psModel.R:logarithmicModel ----
#' @rdname zetaModel
#' @export
logarithmicModel = function() {
  newPsModel(
    model = "logarithmic",
    parameterNames = "pi",
    supportedEngines = c("numerical", "mcmc"),
    subclass = "logarithmicModel"
  )
}

# ---- psModel.R:modelProbabilities.logarithmicModel ----
#' @rdname modelProbabilities
#' @keywords internal
#' @exportS3Method modelProbabilities logarithmicModel
#' @noRd
modelProbabilities.logarithmicModel = function(model, parameters, n, type, ...) {
  logarithmicProbabilities(
    pi = modelParameter(parameters, "pi"),
    n = n,
    type = type
  )
}

# ---- psModel.R:modelLogLikelihood.logarithmicModel ----
#' Evaluate the logarithmic model log likelihood.
#'
#' @param model A `logarithmicModel` descriptor.
#' @param parameters Named parameters containing scalar `pi`.
#' @param data An object of class `psData`.
#' @param ... Additional arguments reserved for future logarithmic controls.
#' @return Scalar log likelihood.
#' @keywords internal
#' @exportS3Method modelLogLikelihood logarithmicModel
#' @noRd
modelLogLikelihood.logarithmicModel = function(model, parameters, data, ...) {
  if (!is(data, "psData")) {
    stop("data must be an object of class psData")
  }

  pi = modelParameter(parameters, "pi")
  validateLogarithmicPi(pi)

  observations = modelObservationData(model, data)
  sum(data$data$rn * dlog(observations, shape = pi, log = TRUE))
}

# ---- fit.R:fitModel.logarithmicModel ----
#' @rdname fitModel
#' @keywords internal
#' @exportS3Method fitModel logarithmicModel
#' @noRd
fitModel.logarithmicModel = function(model, x, ...) {
  result = fitlogDistImpl(x = x, ...)
  result$modelObject = model
  result
}
