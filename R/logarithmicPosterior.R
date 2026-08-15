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
