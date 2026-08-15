#' Fit a one-dimensional numerical posterior for a model.
#'
#' This secondary model dispatch keeps the numerical engine independent of
#' concrete model names. Models opt into the engine by implementing a method
#' that returns a numerical posterior representation.
#'
#' @param model An internal `psModel` object.
#' @param engine A `numericalPosteriorEngine` object.
#' @param x An object of class `psData`.
#' @param prior A prior object created by `makePrior()`.
#' @param ... Model-specific numerical integration controls.
#' @return A `numericalPosteriorRepresentation` object.
#' @keywords internal
#' @noRd
fitNumericalPosteriorModel = function(model, engine, x, prior, ...) {
  UseMethod("fitNumericalPosteriorModel")
}

#' @rdname fitNumericalPosteriorModel
#' @keywords internal
#' @exportS3Method fitNumericalPosteriorModel psModel
#' @noRd
fitNumericalPosteriorModel.psModel = function(model, engine, x, prior, ...) {
  stop(
    "Numerical posterior fitting is not yet implemented for model '",
    model$model,
    "'",
    call. = FALSE
  )
}

#' @rdname fitNumericalPosteriorModel
#' @keywords internal
#' @exportS3Method fitNumericalPosteriorModel zetaModel
#' @noRd
fitNumericalPosteriorModel.zetaModel = function(model, engine, x, prior, ...) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData")
  }
  if (!inherits(prior, "psPrior")) {
    stop("prior must be an object of class psPrior")
  }

  modelObservationData(model, x)
  validateZetaPriorRange(prior$range)

  lower = prior$range[1]
  upper = prior$range[2]

  negLogPosterior = function(shape) {
    -(
      modelLogLikelihood(
        model,
        parameters = list(shape = shape),
        data = x
      ) + prior$logd(shape)
    )
  }

  # The shift by the posterior mode preserves density ratios while reducing
  # underflow when the unnormalised posterior is exponentiated for integrate().
  optimum = optim(
    par = mean(c(lower, upper)),
    fn = negLogPosterior,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper
  )
  logShift = optimum$value

  unnormalisedDensity = Vectorize(function(shape) {
    exp(-negLogPosterior(shape) + logShift)
  })

  normalisingIntegral = integrate(
    unnormalisedDensity,
    lower = lower,
    upper = upper
  )
  logNormalisingConstant = log(normalisingIntegral$value)

  density = Vectorize(function(shape) {
    exp(
      -negLogPosterior(shape) +
        logShift -
        logNormalisingConstant
    )
  })

  meanIntegral = integrate(
    function(shape) {
      shape * density(shape)
    },
    lower = lower,
    upper = upper
  )
  secondMomentIntegral = integrate(
    function(shape) {
      shape^2 * density(shape)
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
      mode = c(shape = unname(optimum$par)),
      normalisingError = normalisingIntegral$abs.error,
      meanError = meanIntegral$abs.error,
      secondMomentError = secondMomentIntegral$abs.error
    )
  )
}

#' Fit through the numerical posterior engine.
#'
#' @param engine A numerical posterior engine.
#' @param model An internal `psModel` object.
#' @param x An object of class `psData`.
#' @param prior A prior object created by `makePrior()`.
#' @param ... Model-specific numerical integration controls.
#' @return A `numericalPosteriorRepresentation` object.
#' @keywords internal
#' @exportS3Method fitPosterior numericalPosteriorEngine
#' @noRd
fitPosterior.numericalPosteriorEngine = function(engine,
                                                   model,
                                                   x,
                                                   prior,
                                                   ...) {
  validateEngineModelPair(engine, model)
  fitNumericalPosteriorModel(
    model = model,
    engine = engine,
    x = x,
    prior = prior,
    ...
  )
}

#' Summarise a numerical posterior representation.
#'
#' @param engine A numerical posterior engine.
#' @param model An internal `psModel` object.
#' @param representation A numerical posterior representation.
#' @param ... Additional arguments reserved for future numerical summaries.
#' @return A data frame of posterior parameter means and standard deviations.
#' @keywords internal
#' @exportS3Method summarisePosterior numericalPosteriorEngine
#' @noRd
summarisePosterior.numericalPosteriorEngine = function(engine,
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
    stop("numerical posterior representation does not match model parameters")
  }

  data.frame(
    parameter = parameterNames,
    estimate = unname(posteriorMean),
    sd = sqrt(pmax(0, diag(posteriorVariance))),
    stringsAsFactors = FALSE
  )
}

#' Extract diagnostics from a numerical posterior representation.
#'
#' @param engine A numerical posterior engine.
#' @param representation A numerical posterior representation.
#' @param ... Additional arguments reserved for future diagnostics.
#' @return A named list of numerical integration diagnostics.
#' @keywords internal
#' @exportS3Method posteriorDiagnostics numericalPosteriorEngine
#' @noRd
posteriorDiagnostics.numericalPosteriorEngine = function(engine,
                                                           representation,
                                                           ...) {
  validatePosteriorRepresentation(representation, engine)
  representation$metadata
}

#' Extract the numerical posterior mean as the fit point estimate.
#'
#' @param engine A numerical posterior engine.
#' @param model An internal `psModel` object.
#' @param representation A numerical posterior representation.
#' @param ... Additional arguments passed to `summarisePosterior()`.
#' @return Named numeric vector of posterior mean parameter estimates.
#' @keywords internal
#' @exportS3Method posteriorPointEstimate numericalPosteriorEngine
#' @noRd
posteriorPointEstimate.numericalPosteriorEngine = function(engine,
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
