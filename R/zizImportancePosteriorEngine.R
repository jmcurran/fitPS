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

#' Fit an importance-sampling posterior for the zero-inflated zeta model.
#'
#' This method preserves the established Laplace-based Gaussian proposal,
#' transformed-coordinate sampling, and normalized importance weights while
#' wrapping the result in the common posterior-engine representation.
#'
#' @param model A `zizModel` descriptor.
#' @param engine An `importancePosteriorEngine` object.
#' @param x An object of class `psData`.
#' @param prior A prior object created by `makePrior()`.
#' @param shape1 First beta-prior shape parameter for `pi`.
#' @param shape2 Second beta-prior shape parameter for `pi`.
#' @param nSamples Number of importance samples.
#' @param proposalScale Positive proposal covariance scale multiplier.
#' @param seed Optional random-number seed.
#' @param start Starting values for `pi` and `shape`.
#' @param laplace Optional precomputed Laplace approximation.
#' @param ... Additional importance controls reserved for future use.
#' @return An `importancePosteriorRepresentation` object.
#' @keywords internal
#' @exportS3Method fitImportancePosteriorModel zizModel
#' @noRd
fitImportancePosteriorModel.zizModel = function(model,
                                                 engine,
                                                 x,
                                                 prior,
                                                 shape1 = 1,
                                                 shape2 = 1,
                                                 nSamples = 5000,
                                                 proposalScale = 2,
                                                 seed = NULL,
                                                 start = c(pi = 0.5, shape = 2),
                                                 laplace = NULL,
                                                 ...) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData")
  }
  if (!inherits(prior, "psPrior")) {
    stop("prior must be an object of class psPrior")
  }

  validateZetaPriorRange(prior$range)

  modelObservationData(model, x)
  approximation = makeZizPosteriorImportance(
    x = x,
    prior = prior,
    shape1 = shape1,
    shape2 = shape2,
    nSamples = nSamples,
    proposalScale = proposalScale,
    seed = seed,
    start = start,
    laplace = laplace
  )

  newPsPosteriorRepresentation(
    engine = engine,
    value = list(
      approximation = approximation,
      mean = approximation$mean,
      variance = approximation$varCov
    ),
    metadata = c(
      list(
        model = model$model,
        shape1 = shape1,
        shape2 = shape2,
        seed = seed
      ),
      approximation$diagnostics
    )
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
