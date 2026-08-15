#' Summarise derived P/S probabilities from a posterior representation.
#'
#' @param engine An internal `psPosteriorEngine` object.
#' @param model An internal `psModel` object.
#' @param representation An internal `psPosteriorRepresentation` object.
#' @param x An object of class `psData`.
#' @param nterms Number of P/S probability terms to summarise.
#' @param level Equal-tailed credible interval level.
#' @param ... Engine-specific controls for derived probability summaries.
#' @return A data frame of posterior P/S probability summaries.
#' @keywords internal
#' @noRd
summarisePosteriorProbabilities = function(engine,
                                            model,
                                            representation,
                                            x,
                                            nterms,
                                            level = 0.95,
                                            ...) {
  UseMethod("summarisePosteriorProbabilities")
}

#' @rdname summarisePosteriorProbabilities
#' @keywords internal
#' @exportS3Method summarisePosteriorProbabilities psPosteriorEngine
#' @noRd
summarisePosteriorProbabilities.psPosteriorEngine = function(engine,
                                                              model,
                                                              representation,
                                                              x,
                                                              nterms,
                                                              level = 0.95,
                                                              ...) {
  validateEngineModelPair(engine, model)
  validatePosteriorRepresentation(representation, engine)
  stop(
    "Posterior probability summaries are not implemented for engine '",
    posteriorEngineName(engine),
    "'",
    call. = FALSE
  )
}

#' Return requested probability indices for a fitted survey type.
#'
#' @param type Survey type, either `"P"` or `"S"`.
#' @param nterms Number of P/S probability terms.
#' @return Integer P/S indices.
#' @keywords internal
#' @noRd
posteriorProbabilityIndices = function(type, nterms) {
  type = normaliseSurveyType(type)
  nterms = as.integer(nterms)

  if (!is.finite(nterms) || length(nterms) != 1L || nterms <= 0L) {
    stop("nterms must be a positive integer")
  }

  if (type == "P") {
    seq.int(0L, nterms - 1L)
  } else {
    seq_len(nterms)
  }
}

#' Summarise posterior probabilities from retained MCMC draws.
#'
#' @inheritParams summarisePosteriorProbabilities
#' @return A data frame of posterior P/S probability summaries.
#' @keywords internal
#' @exportS3Method summarisePosteriorProbabilities mcmcPosteriorEngine
#' @noRd
summarisePosteriorProbabilities.mcmcPosteriorEngine = function(engine,
                                                                model,
                                                                representation,
                                                                x,
                                                                nterms,
                                                                level = 0.95,
                                                                ...) {
  validateEngineModelPair(engine, model)
  validatePosteriorRepresentation(representation, engine)

  chain = representation$value$chain
  parameters = if (is.data.frame(chain) || is.matrix(chain)) {
    as.data.frame(chain)
  } else {
    parameterNames = modelParameterNames(model)
    if (length(parameterNames) != 1L) {
      stop("scalar MCMC chains require a one-parameter model")
    }
    structure(list(chain), names = parameterNames)
  }

  probabilities = modelProbabilities(
    model = model,
    parameters = parameters,
    n = posteriorProbabilityIndices(x$type, nterms),
    type = x$type
  )

  summariseZizProbabilities(
    probabilities = probabilities,
    level = level,
    posteriorMethod = posteriorEngineName(engine)
  )
}

#' Summarise numerical posterior probabilities for a model.
#'
#' @param model An internal `psModel` object.
#' @param engine A numerical posterior engine.
#' @param representation A numerical posterior representation.
#' @param x An object of class `psData`.
#' @param nterms Number of P/S probability terms.
#' @param level Equal-tailed credible interval level.
#' @param ... Additional model-specific controls.
#' @return A data frame of posterior P/S probability summaries.
#' @keywords internal
#' @noRd
summariseNumericalPosteriorProbabilities = function(model,
                                                     engine,
                                                     representation,
                                                     x,
                                                     nterms,
                                                     level = 0.95,
                                                     ...) {
  UseMethod("summariseNumericalPosteriorProbabilities")
}

#' @rdname summariseNumericalPosteriorProbabilities
#' @keywords internal
#' @exportS3Method summariseNumericalPosteriorProbabilities psModel
#' @noRd
summariseNumericalPosteriorProbabilities.psModel = function(model,
                                                             engine,
                                                             representation,
                                                             x,
                                                             nterms,
                                                             level = 0.95,
                                                             ...) {
  stop(
    "Numerical posterior probability summaries are not implemented for model '",
    model$model,
    "'",
    call. = FALSE
  )
}
#' @rdname summarisePosteriorProbabilities
#' @keywords internal
#' @exportS3Method summarisePosteriorProbabilities numericalPosteriorEngine
#' @noRd
summarisePosteriorProbabilities.numericalPosteriorEngine = function(engine,
                                                                     model,
                                                                     representation,
                                                                     x,
                                                                     nterms,
                                                                     level = 0.95,
                                                                     ...) {
  validateEngineModelPair(engine, model)
  validatePosteriorRepresentation(representation, engine)
  summariseNumericalPosteriorProbabilities(
    model = model,
    engine = engine,
    representation = representation,
    x = x,
    nterms = nterms,
    level = level,
    ...
  )
}

#' Summarise posterior probabilities from a Laplace representation.
#'
#' @inheritParams summarisePosteriorProbabilities
#' @param nPosteriorDraws Number of Gaussian working-scale draws used for
#'   derived probability summaries.
#' @param seed Optional random-number seed.
#' @return A data frame of posterior P/S probability summaries.
#' @keywords internal
#' @exportS3Method summarisePosteriorProbabilities laplacePosteriorEngine
#' @noRd
summarisePosteriorProbabilities.laplacePosteriorEngine = function(engine,
                                                                   model,
                                                                   representation,
                                                                   x,
                                                                   nterms,
                                                                   level = 0.95,
                                                                   nPosteriorDraws = 5000,
                                                                   seed = NULL,
                                                                   ...) {
  validateEngineModelPair(engine, model)
  validatePosteriorRepresentation(representation, engine)

  if (!inherits(model, "zizModel")) {
    stop("Laplace probability summaries currently require a zizModel")
  }

  nPosteriorDraws = as.integer(nPosteriorDraws)
  if (!is.finite(nPosteriorDraws) || nPosteriorDraws < 100L) {
    stop("nPosteriorDraws must be at least 100")
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("seed must be NULL or a finite numeric value")
    }
    set.seed(as.integer(seed))
  }

  approximation = representation$value$approximation
  workingDraws = makeZizProposalDraws(
    mean = approximation$modeWorking,
    covariance = approximation$covarianceWorking,
    n = nPosteriorDraws
  )
  thetaDraws = t(apply(workingDraws, 1L, zizWorkingToTheta))
  probabilities = modelProbabilities(
    model = model,
    parameters = as.data.frame(thetaDraws),
    n = posteriorProbabilityIndices(x$type, nterms),
    type = x$type
  )

  summariseZizProbabilities(
    probabilities = probabilities,
    level = level,
    posteriorMethod = posteriorEngineName(engine)
  )
}

#' Summarise posterior probabilities from importance samples.
#'
#' @inheritParams summarisePosteriorProbabilities
#' @return A data frame of posterior P/S probability summaries.
#' @keywords internal
#' @exportS3Method summarisePosteriorProbabilities importancePosteriorEngine
#' @noRd
summarisePosteriorProbabilities.importancePosteriorEngine = function(engine,
                                                                      model,
                                                                      representation,
                                                                      x,
                                                                      nterms,
                                                                      level = 0.95,
                                                                      ...) {
  validateEngineModelPair(engine, model)
  validatePosteriorRepresentation(representation, engine)

  approximation = representation$value$approximation
  samples = approximation$samples
  parameterNames = modelParameterNames(model)
  if (!all(c(parameterNames, "weight") %in% names(samples))) {
    stop("Importance posterior representation is incomplete")
  }

  probabilities = modelProbabilities(
    model = model,
    parameters = samples[parameterNames],
    n = posteriorProbabilityIndices(x$type, nterms),
    type = x$type
  )

  summariseZizProbabilities(
    probabilities = probabilities,
    weights = samples$weight,
    level = level,
    posteriorMethod = posteriorEngineName(engine)
  )
}

#' Return established pre-1.0.7 Bayesian fit fields that must remain compatible.
#'
#' @param model An internal `psModel` object.
#' @param representation An internal posterior representation.
#' @return A named list of compatibility fields.
#' @keywords internal
#' @noRd
establishedBayesianFitFields = function(model, representation) {
  UseMethod("establishedBayesianFitFields")
}

#' @rdname establishedBayesianFitFields
#' @keywords internal
#' @exportS3Method establishedBayesianFitFields psModel
#' @noRd
establishedBayesianFitFields.psModel = function(model, representation) {
  list()
}
#' Finalise a Bayesian fit using the common posterior contract.
#'
#' @param model An internal `psModel` object.
#' @param engine An internal `psPosteriorEngine` object.
#' @param representation An internal posterior representation.
#' @param x An object of class `psData`.
#' @param nterms Number of fitted P/S probability terms.
#' @param level Equal-tailed credible interval level.
#' @param ... Engine-specific controls used by posterior probability summaries.
#' @return A Bayesian object of class `psFit` with a common `psPosterior`.
#' @keywords internal
#' @noRd
finaliseBayesianPsFit = function(model,
                                  engine,
                                  representation,
                                  x,
                                  nterms,
                                  level = 0.95,
                                  ...) {
  validateEngineModelPair(engine, model)
  validatePosteriorRepresentation(representation, engine)

  parameterSummary = summarisePosterior(engine, model, representation)
  pointEstimate = posteriorPointEstimate(engine, model, representation)
  diagnostics = posteriorDiagnostics(engine, representation)

  probabilitySummary = summarisePosteriorProbabilities(
    engine = engine,
    model = model,
    representation = representation,
    x = x,
    nterms = nterms,
    level = level,
    ...
  )

  probabilityIndices = posteriorProbabilityIndices(x$type, nterms)
  fittedMatrix = modelProbabilities(
    model = model,
    parameters = as.list(pointEstimate),
    n = probabilityIndices,
    type = x$type
  )
  fitted = as.numeric(fittedMatrix[1L, ])
  names(fitted) = colnames(fittedMatrix)

  posterior = newPsPosterior(
    method = posteriorEngineName(engine),
    parameters = parameterSummary,
    probabilities = probabilitySummary,
    representation = representation,
    level = level,
    diagnostics = diagnostics,
    model = model$model
  )

  result = c(
    list(
      psData = x,
      fit = list(par = pointEstimate),
      fitted = fitted,
      posterior = posterior,
      model = model$model,
      modelObject = model,
      method = "bayes",
      posteriorMethod = posteriorEngineName(engine)
    ),
    establishedBayesianFitFields(model, representation)
  )

  for (parameterName in names(pointEstimate)) {
    result[[parameterName]] = unname(pointEstimate[[parameterName]])
  }

  class(result) = "psFit"
  result
}

#' Fit a Bayesian model through the shared model/engine architecture.
#'
#' @param model An internal `psModel` object.
#' @param posteriorMethod Posterior approximation method.
#' @param x An object of class `psData`.
#' @param prior A prior object created by `makePrior()`.
#' @param nterms Number of fitted P/S probability terms.
#' @param level Equal-tailed credible interval level.
#' @param ... Engine- and model-specific controls.
#' @return A Bayesian `psFit` object.
#' @keywords internal
#' @noRd
fitBayesianModel = function(model,
                             posteriorMethod,
                             x,
                             prior,
                             nterms,
                             level = 0.95,
                             ...) {
  engine = posteriorEngine(posteriorMethod)
  validateEngineModelPair(engine, model)

  representation = fitPosterior(
    engine = engine,
    model = model,
    x = x,
    prior = prior,
    ...
  )

  finaliseBayesianPsFit(
    model = model,
    engine = engine,
    representation = representation,
    x = x,
    nterms = nterms,
    level = level,
    ...
  )
}
