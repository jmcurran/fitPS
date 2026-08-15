#' Fit an importance-sampling posterior for a model.
#'
#' @param model An internal `psModel` object.
#' @param engine An `importancePosteriorEngine` object.
#' @param x An object of class `psData`.
#' @param prior A prior object created by `makePrior()`.
#' @param ... Model-specific importance-sampling controls.
#' @return An `importancePosteriorRepresentation` object when implemented.
#' @keywords internal
#' @noRd
fitImportancePosteriorModel = function(model, engine, x, prior, ...) {
  UseMethod("fitImportancePosteriorModel")
}

#' @rdname fitImportancePosteriorModel
#' @keywords internal
#' @exportS3Method fitImportancePosteriorModel psModel
#' @noRd
fitImportancePosteriorModel.psModel = function(model, engine, x, prior, ...) {
  stop(
    "importance posterior fitting is not implemented for model '",
    model$model,
    "'",
    call. = FALSE
  )
}
#' Fit through the importance-sampling posterior engine.
#'
#' @param engine An importance posterior engine.
#' @param model An internal `psModel` object.
#' @param x An object of class `psData`.
#' @param prior A prior object created by `makePrior()`.
#' @param ... Model-specific importance controls.
#' @return An `importancePosteriorRepresentation` object.
#' @keywords internal
#' @exportS3Method fitPosterior importancePosteriorEngine
#' @noRd
fitPosterior.importancePosteriorEngine = function(engine,
                                                   model,
                                                   x,
                                                   prior,
                                                   ...) {
  validateEngineModelPair(engine, model)
  fitImportancePosteriorModel(
    model = model,
    engine = engine,
    x = x,
    prior = prior,
    ...
  )
}

#' Summarise an importance posterior representation.
#'
#' @param engine An importance posterior engine.
#' @param model An internal `psModel` object.
#' @param representation An importance posterior representation.
#' @param ... Additional arguments reserved for future summaries.
#' @return A data frame of weighted posterior means and standard deviations.
#' @keywords internal
#' @exportS3Method summarisePosterior importancePosteriorEngine
#' @noRd
summarisePosterior.importancePosteriorEngine = function(engine,
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
    stop("importance posterior representation does not match model parameters")
  }

  data.frame(
    parameter = parameterNames,
    estimate = unname(posteriorMean),
    sd = sqrt(pmax(0, diag(posteriorVariance))),
    stringsAsFactors = FALSE
  )
}

#' Extract diagnostics from an importance posterior representation.
#'
#' @param engine An importance posterior engine.
#' @param representation An importance posterior representation.
#' @param ... Additional arguments reserved for future diagnostics.
#' @return A named list of importance-sampling diagnostics.
#' @keywords internal
#' @exportS3Method posteriorDiagnostics importancePosteriorEngine
#' @noRd
posteriorDiagnostics.importancePosteriorEngine = function(engine,
                                                           representation,
                                                           ...) {
  validatePosteriorRepresentation(representation, engine)
  representation$value$approximation$diagnostics
}

#' Extract the weighted posterior mean as the fit point estimate.
#'
#' @param engine An importance posterior engine.
#' @param model An internal `psModel` object.
#' @param representation An importance posterior representation.
#' @param ... Additional arguments passed to `summarisePosterior()`.
#' @return Named numeric vector of weighted posterior mean estimates.
#' @keywords internal
#' @exportS3Method posteriorPointEstimate importancePosteriorEngine
#' @noRd
posteriorPointEstimate.importancePosteriorEngine = function(engine,
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
