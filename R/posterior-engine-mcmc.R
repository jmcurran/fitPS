#' Fit an MCMC posterior for a model.
#'
#' This secondary model dispatch keeps the MCMC engine independent of concrete
#' model names. Built-in models may preserve specialised samplers, while the
#' `psModel` fallback provides a model-neutral random-walk Metropolis sampler
#' for external models implementing the public Bayesian model contract.
#'
#' @param model A `psModel` object.
#' @param engine An `mcmcPosteriorEngine` object.
#' @param x An object of class `psData`.
#' @param prior A model-specific prior specification.
#' @param nIter Number of retained MCMC iterations.
#' @param nBurnIn Number of burn-in iterations.
#' @param proposalScale Positive random-walk proposal standard deviation on the
#'   unconstrained scale. Supply one value or a named value for every model parameter.
#' @param seed Optional integer random seed.
#' @param ... Model-specific Bayesian controls passed to the public model
#'   contract.
#' @return An `mcmcPosteriorRepresentation` object.
#' @keywords internal
#' @noRd
fitMcmcPosteriorModel = function(model, engine, x, prior, ...) {
  UseMethod("fitMcmcPosteriorModel")
}

validateGenericMcmcIterations = function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value <= 0 || value != floor(value)) {
    stop(name, " must be one positive integer", call. = FALSE)
  }
  as.integer(value)
}

normaliseGenericProposalScale = function(model, proposalScale) {
  parameterNames = modelParameterNames(model)

  if (!is.numeric(proposalScale) || any(!is.finite(proposalScale)) ||
      any(proposalScale <= 0)) {
    stop("proposalScale must contain positive finite numeric values", call. = FALSE)
  }

  if (length(proposalScale) == 1L) {
    proposalScale = rep(proposalScale, length(parameterNames))
    names(proposalScale) = parameterNames
    return(proposalScale)
  }

  if (length(proposalScale) != length(parameterNames) ||
      is.null(names(proposalScale)) ||
      !setequal(names(proposalScale), parameterNames)) {
    stop(
      "proposalScale must be one value or a named numeric vector matching model parameters",
      call. = FALSE
    )
  }

  proposalScale[parameterNames]
}

#' @rdname fitMcmcPosteriorModel
#' @keywords internal
#' @exportS3Method fitMcmcPosteriorModel psModel
#' @noRd
#' @importFrom stats cov rnorm runif
fitMcmcPosteriorModel.psModel = function(model,
                                         engine,
                                         x,
                                         prior,
                                         nIter = 2000L,
                                         nBurnIn = 500L,
                                         proposalScale = 0.25,
                                         seed = NULL,
                                         ...) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData", call. = FALSE)
  }

  nIter = validateGenericMcmcIterations(nIter, "nIter")
  nBurnIn = validateGenericMcmcIterations(nBurnIn, "nBurnIn")
  proposalScale = normaliseGenericProposalScale(model, proposalScale)

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("seed must be NULL or one finite numeric value", call. = FALSE)
    }
    set.seed(as.integer(seed))
  }

  controls = modelBayesControl(
    model = model,
    x = x,
    engine = engine,
    prior = prior,
    ...
  )
  if (!is.list(controls) || is.null(controls$start)) {
    stop("modelBayesControl() must return a list containing start", call. = FALSE)
  }

  currentNatural = validateModelParameterVector(model, controls$start, "start")
  currentUnconstrained = modelToUnconstrained(model, currentNatural)
  currentUnconstrained = validateModelParameterVector(model, currentUnconstrained, "unconstrained start")
  parameterNames = modelParameterNames(model)

  logUnconstrainedPosterior = function(unconstrained) {
    names(unconstrained) = parameterNames

    natural = tryCatch(
      modelFromUnconstrained(model, unconstrained),
      error = function(error) NULL
    )
    if (is.null(natural)) {
      return(-Inf)
    }

    natural = tryCatch(
      validateModelParameterVector(model, natural, "natural parameters"),
      error = function(error) NULL
    )
    if (is.null(natural)) {
      return(-Inf)
    }

    logLikelihood = modelLogLikelihood(
      model = model,
      parameters = natural,
      data = x
    )
    logPrior = modelLogPrior(
      model = model,
      parameters = natural,
      prior = prior,
      ...
    )
    logJacobian = modelLogJacobian(model, unconstrained)

    values = c(logLikelihood, logPrior, logJacobian)
    if (!is.numeric(values) || length(values) != 3L || any(is.na(values))) {
      return(-Inf)
    }
    if (any(values == -Inf)) {
      return(-Inf)
    }
    if (any(!is.finite(values))) {
      return(-Inf)
    }

    sum(values)
  }

  currentLogPosterior = logUnconstrainedPosterior(currentUnconstrained)
  if (!is.finite(currentLogPosterior)) {
    stop("posterior is not finite at the model-supplied starting values", call. = FALSE)
  }

  nTotal = nIter + nBurnIn
  chain = matrix(
    NA_real_,
    nrow = nIter,
    ncol = length(parameterNames),
    dimnames = list(NULL, parameterNames)
  )
  accepted = 0L

  for (iteration in seq_len(nTotal)) {
    proposedUnconstrained = currentUnconstrained + rnorm(
      length(parameterNames),
      mean = 0,
      sd = proposalScale
    )
    names(proposedUnconstrained) = parameterNames
    proposedLogPosterior = logUnconstrainedPosterior(proposedUnconstrained)

    if (is.finite(proposedLogPosterior) &&
        (proposedLogPosterior >= currentLogPosterior ||
          log(runif(1L)) < proposedLogPosterior - currentLogPosterior)) {
      currentUnconstrained = proposedUnconstrained
      currentLogPosterior = proposedLogPosterior
      accepted = accepted + 1L
    }

    if (iteration > nBurnIn) {
      chain[iteration - nBurnIn, ] = modelFromUnconstrained(model, currentUnconstrained)
    }
  }

  chain = as.data.frame(chain)
  posteriorMean = vapply(chain, mean, numeric(1L))
  posteriorVariance = if (length(parameterNames) == 1L) {
    matrix(
      var(chain[[1L]]),
      nrow = 1L,
      dimnames = list(parameterNames, parameterNames)
    )
  } else {
    cov(chain)
  }

  newPsPosteriorRepresentation(
    engine = engine,
    value = list(
      chain = chain,
      mean = posteriorMean,
      variance = posteriorVariance
    ),
    metadata = list(
      model = model$model,
      nIter = nIter,
      nBurnIn = nBurnIn,
      proposalScale = proposalScale,
      acceptance = accepted / nTotal,
      seed = seed,
      generic = TRUE
    )
  )
}

#' Fit through the MCMC posterior engine.
#'
#' @param engine An MCMC posterior engine.
#' @param model A `psModel` object.
#' @param x An object of class `psData`.
#' @param prior A model-specific prior specification.
#' @param ... Model-specific MCMC controls.
#' @return An `mcmcPosteriorRepresentation` object.
#' @keywords internal
#' @exportS3Method fitPosterior mcmcPosteriorEngine
#' @noRd
fitPosterior.mcmcPosteriorEngine = function(engine, model, x, prior, ...) {
  validateEngineModelPair(engine, model)
  fitMcmcPosteriorModel(
    model = model,
    engine = engine,
    x = x,
    prior = prior,
    ...
  )
}
#' Summarise an MCMC posterior representation.
#'
#' @param engine An MCMC posterior engine.
#' @param model An internal `psModel` object.
#' @param representation An MCMC posterior representation.
#' @param ... Additional arguments reserved for future MCMC summaries.
#' @return A data frame of posterior parameter means and standard deviations.
#' @keywords internal
#' @exportS3Method summarisePosterior mcmcPosteriorEngine
#' @noRd
summarisePosterior.mcmcPosteriorEngine = function(engine,
                                                   model,
                                                   representation,
                                                   ...) {
  validateEngineModelPair(engine, model)
  validatePosteriorRepresentation(representation, engine)

  parameterNames = modelParameterNames(model)
  posteriorMean = representation$value$mean
  posteriorVariance = representation$value$variance

  varianceDimensions = dim(posteriorVariance)
  expectedDimensions = c(length(parameterNames), length(parameterNames))
  if (length(posteriorMean) != length(parameterNames) ||
      is.null(varianceDimensions) ||
      !identical(as.integer(varianceDimensions), as.integer(expectedDimensions))) {
    stop("MCMC posterior representation does not match model parameters")
  }

  data.frame(
    parameter = parameterNames,
    estimate = unname(posteriorMean),
    sd = sqrt(pmax(0, diag(posteriorVariance))),
    stringsAsFactors = FALSE
  )
}

#' Extract diagnostics from an MCMC posterior representation.
#'
#' @param engine An MCMC posterior engine.
#' @param representation An MCMC posterior representation.
#' @param ... Additional arguments reserved for future diagnostics.
#' @return A named list of MCMC diagnostics and run settings.
#' @keywords internal
#' @exportS3Method posteriorDiagnostics mcmcPosteriorEngine
#' @noRd
posteriorDiagnostics.mcmcPosteriorEngine = function(engine,
                                                     representation,
                                                     ...) {
  validatePosteriorRepresentation(representation, engine)
  representation$metadata
}

#' Extract the MCMC posterior mean as the fit point estimate.
#'
#' @param engine An MCMC posterior engine.
#' @param model An internal `psModel` object.
#' @param representation An MCMC posterior representation.
#' @param ... Additional arguments passed to `summarisePosterior()`.
#' @return Named numeric vector of posterior mean parameter estimates.
#' @keywords internal
#' @exportS3Method posteriorPointEstimate mcmcPosteriorEngine
#' @noRd
posteriorPointEstimate.mcmcPosteriorEngine = function(engine,
                                                       model,
                                                       representation,
                                                       ...) {
  summary = summarisePosterior(
    engine = engine,
    model = model,
    representation = representation,
    ...
  )
  result = summary$estimate
  names(result) = summary$parameter
  result
}
