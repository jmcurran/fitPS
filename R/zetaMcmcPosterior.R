#' Fit an MCMC posterior for a model.
#'
#' This secondary model dispatch keeps the MCMC engine independent of concrete
#' model names. Models opt into the engine by implementing a method that
#' returns an MCMC posterior representation.
#'
#' @param model An internal `psModel` object.
#' @param engine An `mcmcPosteriorEngine` object.
#' @param x An object of class `psData`.
#' @param prior A prior object created by `makePrior()`.
#' @param ... Model-specific MCMC controls.
#' @return An `mcmcPosteriorRepresentation` object.
#' @keywords internal
#' @noRd
fitMcmcPosteriorModel = function(model, engine, x, prior, ...) {
  UseMethod("fitMcmcPosteriorModel")
}

#' @rdname fitMcmcPosteriorModel
#' @keywords internal
#' @exportS3Method fitMcmcPosteriorModel psModel
#' @noRd
fitMcmcPosteriorModel.psModel = function(model, engine, x, prior, ...) {
  stop(
    "MCMC posterior fitting is not yet implemented for model '",
    model$model,
    "'",
    call. = FALSE
  )
}

#' @rdname fitMcmcPosteriorModel
#' @keywords internal
#' @exportS3Method fitMcmcPosteriorModel zetaModel
#' @noRd
fitMcmcPosteriorModel.zetaModel = function(model,
                                            engine,
                                            x,
                                            prior,
                                            shape0 = 2,
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

  validateZetaShape(shape0, "shape0")
  validatePriorRange(prior$range)
  modelObservationData(model, x)

  if (nIter < 1000) {
    warning(
      "The number of samples from the MCMC chain really should be 1000 or higher."
    )
  }
  if (nIter <= 0 || nBurnIn <= 0) {
    stop("nIter and nBurnIn must be greater than zero.")
  }

  lower = prior$range[1]
  upper = prior$range[2]
  nTotal = nIter + nBurnIn

  # Draw the complete proposal and acceptance streams before the loop. Keeping
  # this ordering preserves the legacy RNG sequence for a fixed set.seed().
  proposals = runif(nTotal, lower, upper)
  logUniforms = log(runif(nTotal))
  chain = numeric(nIter)

  logPosterior = function(shape) {
    modelLogLikelihood(
      model = model,
      parameters = list(shape = shape),
      data = x
    ) + prior$logd(shape)
  }

  currentLogPosterior = logPosterior(shape0)
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

  i = 1
  while (i <= nTotal) {
    proposedShape = proposals[i]
    proposedLogPosterior = logPosterior(proposedShape)

    if (proposedLogPosterior > currentLogPosterior ||
        logUniforms[i] < (proposedLogPosterior - currentLogPosterior)) {
      shape0 = proposedShape
      currentLogPosterior = proposedLogPosterior
    }

    if (i > nBurnIn) {
      chain[i - nBurnIn] = shape0
    }

    i = i + 1
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
      mean = c(shape = posteriorMean),
      variance = matrix(
        posteriorVariance,
        nrow = 1L,
        dimnames = list("shape", "shape")
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

#' Fit through the MCMC posterior engine.
#'
#' @param engine An MCMC posterior engine.
#' @param model An internal `psModel` object.
#' @param x An object of class `psData`.
#' @param prior A prior object created by `makePrior()`.
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
