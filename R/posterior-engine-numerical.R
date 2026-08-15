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
#' @importFrom stats setNames
#' @noRd
fitNumericalPosteriorModel.psModel = function(model, engine, x, prior, ...) {
  parameterNames = modelParameterNames(model)
  if (length(parameterNames) != 1L) {
    stop(
      "generic numerical posterior fitting currently requires a one-parameter model",
      call. = FALSE
    )
  }

  parameterName = parameterNames[[1L]]
  control = modelBayesControl(model, x, engine, prior, ...)
  requiredControls = c("start", "lower", "upper")
  if (!is.list(control) || !all(requiredControls %in% names(control))) {
    stop(
      "generic numerical fitting requires modelBayesControl() to return start, lower, and upper",
      call. = FALSE
    )
  }

  validateScalarControl = function(value, name, allowInfinite = FALSE) {
    if (!is.numeric(value) || length(value) != 1L ||
        is.null(names(value)) || !identical(names(value), parameterName) ||
        (!allowInfinite && !is.finite(value))) {
      stop(name, " must be one named numeric value for ", parameterName, call. = FALSE)
    }
    if (allowInfinite && is.na(value)) {
      stop(name, " must not be NA", call. = FALSE)
    }
    unname(value)
  }

  initial = validateScalarControl(control$start, "start")
  lower = validateScalarControl(control$lower, "lower", allowInfinite = TRUE)
  upper = validateScalarControl(control$upper, "upper", allowInfinite = TRUE)
  if (lower >= upper || initial <= lower || initial >= upper) {
    stop("numerical start and bounds must satisfy lower < start < upper", call. = FALSE)
  }

  logPosterior = function(parameterValue) {
    parameters = parameterValue
    names(parameters) = parameterName
    logPrior = modelLogPrior(model, parameters = parameters, prior = prior, ...)
    if (!is.numeric(logPrior) || length(logPrior) != 1L || is.na(logPrior)) {
      stop("modelLogPrior() must return one numeric value", call. = FALSE)
    }
    if (!is.finite(logPrior)) {
      return(-Inf)
    }

    logLikelihood = modelLogLikelihood(model, parameters = parameters, data = x)
    if (!is.numeric(logLikelihood) || length(logLikelihood) != 1L || is.na(logLikelihood)) {
      stop("modelLogLikelihood() must return one numeric value", call. = FALSE)
    }
    if (!is.finite(logLikelihood)) {
      return(-Inf)
    }
    logLikelihood + logPrior
  }

  negLogPosterior = function(parameterValue) {
    value = logPosterior(parameterValue)
    if (!is.finite(value)) {
      return(.Machine$double.xmax^0.5)
    }
    -value
  }

  optimisationLower = lower
  optimisationUpper = upper
  if (is.finite(lower) && is.finite(upper)) {
    boundaryInset = sqrt(.Machine$double.eps) * (upper - lower)
    optimisationLower = lower + boundaryInset
    optimisationUpper = upper - boundaryInset
  }

  optimisationStart = if (is.finite(lower) && is.finite(upper)) {
    mean(c(lower, upper))
  } else {
    initial
  }

  optimum = optim(
    par = optimisationStart,
    fn = negLogPosterior,
    method = "L-BFGS-B",
    lower = optimisationLower,
    upper = optimisationUpper
  )
  if (!identical(optimum$convergence, 0L)) {
    stop("generic numerical posterior optimisation failed", call. = FALSE)
  }
  logShift = optimum$value

  unnormalisedDensity = Vectorize(function(parameterValue) {
    value = logPosterior(parameterValue)
    if (!is.finite(value)) {
      return(0)
    }
    exp(value + logShift)
  })
  normalisingIntegral = integrate(
    unnormalisedDensity,
    lower = lower,
    upper = upper
  )
  logNormalisingConstant = log(normalisingIntegral$value)

  density = Vectorize(function(parameterValue) {
    value = logPosterior(parameterValue)
    if (!is.finite(value)) {
      return(0)
    }
    exp(value + logShift - logNormalisingConstant)
  })
  meanIntegral = integrate(
    function(parameterValue) {
      parameterValue * density(parameterValue)
    },
    lower = lower,
    upper = upper
  )
  secondMomentIntegral = integrate(
    function(parameterValue) {
      parameterValue^2 * density(parameterValue)
    },
    lower = lower,
    upper = upper
  )

  posteriorMean = meanIntegral$value
  posteriorVariance = secondMomentIntegral$value - posteriorMean^2
  meanValue = posteriorMean
  names(meanValue) = parameterName
  variance = matrix(
    posteriorVariance,
    nrow = 1L,
    dimnames = list(parameterName, parameterName)
  )

  newPsPosteriorRepresentation(
    engine = engine,
    value = list(
      density = density,
      mean = meanValue,
      variance = variance
    ),
    metadata = list(
      model = model$model,
      bounds = c(lower = lower, upper = upper),
      mode = setNames(unname(optimum$par), parameterName),
      normalisingError = normalisingIntegral$abs.error,
      meanError = meanIntegral$abs.error,
      secondMomentError = secondMomentIntegral$abs.error,
      generic = TRUE
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
