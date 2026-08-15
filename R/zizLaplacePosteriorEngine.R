#' Fit a Laplace posterior for a model.
#'
#' @param model An internal `psModel` object.
#' @param engine A `laplacePosteriorEngine` object.
#' @param x An object of class `psData`.
#' @param prior A prior object created by `makePrior()`.
#' @param ... Model-specific Laplace controls.
#' @return A `laplacePosteriorRepresentation` object when implemented.
#' @keywords internal
#' @noRd
fitLaplacePosteriorModel = function(model, engine, x, prior, ...) {
  UseMethod("fitLaplacePosteriorModel")
}

#' @rdname fitLaplacePosteriorModel
#' @keywords internal
#' @exportS3Method fitLaplacePosteriorModel psModel
#' @noRd
fitLaplacePosteriorModel.psModel = function(model, engine, x, prior, ...) {
  stop(
    "Laplace posterior fitting is not implemented for model '",
    model$model,
    "'",
    call. = FALSE
  )
}

#' Fit a Laplace posterior for the zero-inflated zeta model.
#'
#' This method preserves the established transformed-coordinate optimisation,
#' Hessian calculation, and delta-method covariance while wrapping the result
#' in the common posterior-engine representation.
#'
#' @param model A `zizModel` descriptor.
#' @param engine A `laplacePosteriorEngine` object.
#' @param x An object of class `psData`.
#' @param prior A prior object created by `makePrior()`.
#' @param shape1 First beta-prior shape parameter for `pi`.
#' @param shape2 Second beta-prior shape parameter for `pi`.
#' @param start Starting values for `pi` and `shape`.
#' @param ... Additional Laplace controls reserved for future use.
#' @return A `laplacePosteriorRepresentation` object.
#' @keywords internal
#' @exportS3Method fitLaplacePosteriorModel zizModel
#' @noRd
fitLaplacePosteriorModel.zizModel = function(model,
                                              engine,
                                              x,
                                              prior,
                                              shape1 = 1,
                                              shape2 = 1,
                                              start = c(pi = 0.5, shape = 2),
                                              ...) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData")
  }
  if (!inherits(prior, "psPrior")) {
    stop("prior must be an object of class psPrior")
  }

  validateZetaPriorRange(prior$range)

  modelObservationData(model, x)
  approximation = makeZizPosteriorLaplace(
    x = x,
    prior = prior,
    shape1 = shape1,
    shape2 = shape2,
    start = start
  )

  newPsPosteriorRepresentation(
    engine = engine,
    value = list(
      approximation = approximation,
      mean = approximation$mode,
      variance = approximation$varCov
    ),
    metadata = list(
      model = model$model,
      convergence = approximation$convergence,
      logPosteriorMode = approximation$logPosteriorMode,
      shape1 = shape1,
      shape2 = shape2
    )
  )
}

#' Fit through the Laplace posterior engine.
#'
#' @param engine A Laplace posterior engine.
#' @param model An internal `psModel` object.
#' @param x An object of class `psData`.
#' @param prior A prior object created by `makePrior()`.
#' @param ... Model-specific Laplace controls.
#' @return A `laplacePosteriorRepresentation` object.
#' @keywords internal
#' @exportS3Method fitPosterior laplacePosteriorEngine
#' @noRd
fitPosterior.laplacePosteriorEngine = function(engine,
                                                model,
                                                x,
                                                prior,
                                                ...) {
  validateEngineModelPair(engine, model)
  fitLaplacePosteriorModel(
    model = model,
    engine = engine,
    x = x,
    prior = prior,
    ...
  )
}

#' Summarise a Laplace posterior representation.
#'
#' @param engine A Laplace posterior engine.
#' @param model An internal `psModel` object.
#' @param representation A Laplace posterior representation.
#' @param ... Additional arguments reserved for future summaries.
#' @return A data frame of posterior mode estimates and local standard deviations.
#' @keywords internal
#' @exportS3Method summarisePosterior laplacePosteriorEngine
#' @noRd
summarisePosterior.laplacePosteriorEngine = function(engine,
                                                      model,
                                                      representation,
                                                      ...) {
  validateEngineModelPair(engine, model)
  validatePosteriorRepresentation(representation, engine)

  parameterNames = modelParameterNames(model)
  posteriorMode = representation$value$mean
  posteriorVariance = representation$value$variance
  varianceDimensions = dim(posteriorVariance)
  expectedDimensions = c(length(parameterNames), length(parameterNames))

  if (length(posteriorMode) != length(parameterNames) ||
      is.null(varianceDimensions) ||
      !identical(as.integer(varianceDimensions), as.integer(expectedDimensions))) {
    stop("Laplace posterior representation does not match model parameters")
  }

  data.frame(
    parameter = parameterNames,
    estimate = unname(posteriorMode),
    sd = sqrt(pmax(0, diag(posteriorVariance))),
    stringsAsFactors = FALSE
  )
}

#' Extract diagnostics from a Laplace posterior representation.
#'
#' @param engine A Laplace posterior engine.
#' @param representation A Laplace posterior representation.
#' @param ... Additional arguments reserved for future diagnostics.
#' @return A named list of Laplace optimisation diagnostics.
#' @keywords internal
#' @exportS3Method posteriorDiagnostics laplacePosteriorEngine
#' @noRd
posteriorDiagnostics.laplacePosteriorEngine = function(engine,
                                                        representation,
                                                        ...) {
  validatePosteriorRepresentation(representation, engine)
  representation$metadata
}

#' Extract the Laplace posterior mode as the fit point estimate.
#'
#' @param engine A Laplace posterior engine.
#' @param model An internal `psModel` object.
#' @param representation A Laplace posterior representation.
#' @param ... Additional arguments passed to `summarisePosterior()`.
#' @return Named numeric vector of posterior mode estimates.
#' @keywords internal
#' @exportS3Method posteriorPointEstimate laplacePosteriorEngine
#' @noRd
posteriorPointEstimate.laplacePosteriorEngine = function(engine,
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
