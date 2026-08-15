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
#' @importFrom cubature hcubature
#' @importFrom stats setNames
#' @noRd
fitNumericalPosteriorModel.psModel = function(model, engine, x, prior, ...) {
  parameterNames = modelParameterNames(model)
  dimension = length(parameterNames)

  if (dimension == 1L) {
    return(fitNumericalPosterior1d(model, engine, x, prior, ...))
  }
  if (dimension == 2L) {
    return(fitNumericalPosterior2d(model, engine, x, prior, ...))
  }

  stop(
    "Numerical posterior fitting in fitPS supports models with at most two parameters. ",
    "Use the MCMC posterior engine for models with three or more parameters.",
    call. = FALSE
  )
}

#' Fit a one-dimensional numerical posterior.
#'
#' @keywords internal
#' @noRd
fitNumericalPosterior1d = function(model, engine, x, prior, ...) {
  parameterNames = modelParameterNames(model)
  parameterName = parameterNames[[1L]]
  control = modelBayesControl(model, x, engine, prior, ...)
  requiredControls = c("start", "lower", "upper")
  if (!is.list(control) || !all(requiredControls %in% names(control))) {
    stop(
      "generic one-dimensional numerical fitting requires modelBayesControl() to return start, lower, and upper",
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
      dimension = 1L,
      integrationMethod = "integrate",
      bounds = c(lower = lower, upper = upper),
      mode = setNames(unname(optimum$par), parameterName),
      normalisingError = normalisingIntegral$abs.error,
      meanError = meanIntegral$abs.error,
      secondMomentError = secondMomentIntegral$abs.error,
      generic = TRUE
    )
  )
}

#' Fit a two-dimensional numerical posterior by adaptive cubature.
#'
#' The model supplies finite natural-scale lower and upper bounds through
#' `modelBayesControl()`. The engine integrates the posterior directly over
#' that rectangle, avoiding an additional transformation that can make a
#' concentrated posterior difficult for adaptive cubature to locate.
#'
#' @keywords internal
#' @noRd
fitNumericalPosterior2d = function(model,
                                    engine,
                                    x,
                                    prior,
                                    tol = 1e-5,
                                    maxEval = 0L,
                                    summaryGridSize = 41L,
                                    ...) {
  parameterNames = modelParameterNames(model)
  control = modelBayesControl(model, x, engine, prior, ...)
  requiredControls = c("start", "lower", "upper")
  if (!is.list(control) || !all(requiredControls %in% names(control))) {
    stop(
      "generic two-dimensional numerical fitting requires modelBayesControl() to return start, lower, and upper",
      call. = FALSE
    )
  }

  validateVectorControl = function(value, name) {
    if (!is.numeric(value) || length(value) != 2L || is.null(names(value)) ||
        !setequal(names(value), parameterNames) || any(!is.finite(value))) {
      stop(name, " must contain two finite named model parameters", call. = FALSE)
    }
    value[parameterNames]
  }

  initial = validateVectorControl(control$start, "start")
  lower = validateVectorControl(control$lower, "lower")
  upper = validateVectorControl(control$upper, "upper")
  if (any(lower >= upper) || any(initial <= lower) || any(initial >= upper)) {
    stop(
      "two-dimensional numerical start and bounds must satisfy lower < start < upper for every parameter",
      call. = FALSE
    )
  }

  tol = as.numeric(tol)
  maxEval = as.integer(maxEval)
  summaryGridSize = as.integer(summaryGridSize)
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("tol must be one positive finite number", call. = FALSE)
  }
  if (length(maxEval) != 1L || is.na(maxEval) || maxEval < 0L) {
    stop("maxEval must be a non-negative integer", call. = FALSE)
  }
  if (length(summaryGridSize) != 1L || is.na(summaryGridSize) || summaryGridSize < 11L) {
    stop("summaryGridSize must be at least 11", call. = FALSE)
  }

  naturalLogPosterior = function(parameters) {
    parameters = as.numeric(parameters)
    names(parameters) = parameterNames
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

  boundaryInset = sqrt(.Machine$double.eps) * pmax(1, upper - lower)
  optimisationLower = lower + boundaryInset
  optimisationUpper = upper - boundaryInset
  optimisationStart = pmin(
    optimisationUpper,
    pmax(optimisationLower, initial)
  )
  optimum = optim(
    par = optimisationStart,
    fn = function(parameters) {
      value = naturalLogPosterior(parameters)
      if (!is.finite(value)) {
        return(.Machine$double.xmax^0.5)
      }
      -value
    },
    method = "L-BFGS-B",
    lower = optimisationLower,
    upper = optimisationUpper
  )
  if (!identical(optimum$convergence, 0L)) {
    stop("generic two-dimensional numerical posterior optimisation failed", call. = FALSE)
  }
  mode = setNames(unname(optimum$par), parameterNames)
  logScale = -optimum$value

  unitToNatural = function(unitParameters) {
    unitParameters = as.numeric(unitParameters)
    belowMode = unitParameters <= 0.5
    parameters = numeric(2L)
    jacobian = numeric(2L)

    parameters[belowMode] = lower[belowMode] +
      2 * unitParameters[belowMode] * (mode[belowMode] - lower[belowMode])
    jacobian[belowMode] = 2 * (mode[belowMode] - lower[belowMode])

    aboveMode = !belowMode
    parameters[aboveMode] = mode[aboveMode] +
      2 * (unitParameters[aboveMode] - 0.5) *
        (upper[aboveMode] - mode[aboveMode])
    jacobian[aboveMode] = 2 * (upper[aboveMode] - mode[aboveMode])

    names(parameters) = parameterNames
    list(parameters = parameters, jacobian = prod(jacobian))
  }

  momentIntegrand = function(unitParameters) {
    transformed = unitToNatural(unitParameters)
    parameters = transformed$parameters
    logValue = naturalLogPosterior(parameters)
    if (!is.finite(logValue) || !is.finite(transformed$jacobian) ||
        transformed$jacobian <= 0) {
      return(rep(0, 7L))
    }
    weight = exp(logValue - logScale) * transformed$jacobian
    c(
      weight,
      parameters * weight,
      parameters[[1L]]^2 * weight,
      parameters[[1L]] * parameters[[2L]] * weight,
      parameters[[2L]]^2 * weight,
      modelDeviance(model, parameters, x) * weight
    )
  }

  cubatureResult = suppressWarnings(hcubature(
    f = momentIntegrand,
    lowerLimit = c(0, 0),
    upperLimit = c(1, 1),
    fDim = 7L,
    tol = tol,
    maxEval = maxEval
  ))
  if (!identical(cubatureResult$returnCode, 0L)) {
    stop(
      "two-dimensional numerical posterior cubature did not converge; use MCMC, revise numerical bounds, or relax numerical controls",
      call. = FALSE
    )
  }

  integrals = cubatureResult$integral
  normalisingScaled = integrals[[1L]]
  if (!is.finite(normalisingScaled) || normalisingScaled <= 0) {
    stop(
      "two-dimensional numerical posterior could not be normalized within the model-supplied bounds",
      call. = FALSE
    )
  }
  posteriorMean = integrals[2:3] / normalisingScaled
  names(posteriorMean) = parameterNames
  secondMoments = matrix(
    c(integrals[[4L]], integrals[[5L]], integrals[[5L]], integrals[[6L]]),
    nrow = 2L
  ) / normalisingScaled
  posteriorVariance = secondMoments - tcrossprod(posteriorMean)
  dimnames(posteriorVariance) = list(parameterNames, parameterNames)

  logNormalisingConstant = log(normalisingScaled) + logScale
  density = function(parameters) {
    if (is.list(parameters)) {
      parameters = unlist(parameters, use.names = TRUE)
    }
    if (!is.numeric(parameters) || length(parameters) != 2L) {
      stop("two-dimensional numerical density requires two model parameters", call. = FALSE)
    }
    if (is.null(names(parameters))) {
      names(parameters) = parameterNames
    }
    parameters = parameters[parameterNames]
    if (any(!is.finite(parameters)) || any(parameters < lower) || any(parameters > upper)) {
      return(0)
    }
    value = naturalLogPosterior(parameters)
    if (!is.finite(value)) {
      return(0)
    }
    exp(value - logNormalisingConstant)
  }

  gridAxes = lapply(seq_along(parameterNames), function(index) {
    seq(lower[[index]], upper[[index]], length.out = summaryGridSize)
  })
  parameterGrid = expand.grid(gridAxes[[1L]], gridAxes[[2L]])
  names(parameterGrid) = parameterNames
  logWeights = apply(parameterGrid, 1L, naturalLogPosterior)
  finiteWeights = is.finite(logWeights)
  if (!any(finiteWeights)) {
    stop("Unable to construct a two-dimensional posterior summary grid", call. = FALSE)
  }
  gridShift = max(logWeights[finiteWeights])
  weights = exp(logWeights - gridShift)
  weights[!is.finite(weights)] = 0
  weightTotal = sum(weights)
  if (!is.finite(weightTotal) || weightTotal <= 0) {
    stop("Unable to normalize the two-dimensional posterior summary grid", call. = FALSE)
  }
  weights = weights / weightTotal

  newPsPosteriorRepresentation(
    engine = engine,
    value = list(
      density = density,
      mean = posteriorMean,
      variance = posteriorVariance,
      grid = list(
        parameters = parameterGrid,
        weights = weights
      )
    ),
    metadata = list(
      model = model$model,
      dimension = 2L,
      integrationMethod = "hcubature",
      bounds = list(lower = lower, upper = upper),
      mode = mode,
      cubatureError = cubatureResult$error,
      cubatureReturnCode = cubatureResult$returnCode,
      expectedDeviance = integrals[[7L]] / normalisingScaled,
      functionEvaluations = cubatureResult$functionEvaluations,
      tolerance = tol,
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
