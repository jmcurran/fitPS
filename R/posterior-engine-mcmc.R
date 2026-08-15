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
