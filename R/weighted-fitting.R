#' Validate weights for internal weighted fitting
#'
#' @param x An object of class `psData`.
#' @param weights Positive finite fitting weights aligned with occupied support
#'   rows in `x`.
#' @return The validated numeric weight vector.
#' @keywords internal
#' @noRd
validateFittingWeights = function(x, weights) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData")
  }
  if (!is.numeric(weights) || length(weights) != nrow(x$data)) {
    stop("weights must be numeric and aligned with the occupied support rows")
  }
  if (any(!is.finite(weights)) || any(weights <= 0)) {
    stop("weights must contain only finite positive values")
  }

  as.numeric(weights)
}

#' Construct an internal weighted survey copy
#'
#' This helper deliberately bypasses [makePSData()] because its `rn` values are
#' inferential fitting weights rather than observed integer survey counts. The
#' returned object is used only inside weighted likelihood calculations.
#'
#' @param x An object of class `psData`.
#' @param weights Validated or candidate fitting weights.
#' @return A copy of `x` whose `rn` column contains the supplied real weights.
#' @keywords internal
#' @noRd
weightedSurveyData = function(x, weights) {
  weights = validateFittingWeights(x, weights)
  result = x
  result$data$rn = weights
  result
}

#' Return optimisation controls for an internal built-in weighted MLE
#'
#' @param model A built-in fitPS model descriptor.
#' @param ... Optional model-specific starting values.
#' @return A list containing `start`, `lower`, and `upper` vectors.
#' @keywords internal
#' @noRd
weightedMleControl = function(model, ...) {
  dots = list(...)

  if (inherits(model, "zetaModel")) {
    start = if ("start" %in% names(dots)) {
      dots$start
    } else if ("shape" %in% names(dots)) {
      dots$shape
    } else {
      2
    }
    validateZetaShape(start, "start")
    return(list(
      start = c(shape = start),
      lower = c(shape = 1 + sqrt(.Machine$double.eps)),
      upper = c(shape = Inf)
    ))
  }

  if (inherits(model, "zizModel")) {
    start = if ("start" %in% names(dots)) {
      dots$start
    } else {
      c(pi = 0.5, shape = 2)
    }
    if (!is.numeric(start) || length(start) != 2L || any(!is.finite(start))) {
      stop("start must contain finite pi and shape values")
    }
    start = as.numeric(start)
    if (start[1L] <= 0 || start[1L] >= 1) {
      stop("The starting value for pi must be in (0, 1)")
    }
    validateZetaShape(start[2L], "start shape")
    return(list(
      start = c(pi = start[1L], shape = start[2L]),
      lower = c(
        pi = sqrt(.Machine$double.eps),
        shape = 1 + sqrt(.Machine$double.eps)
      ),
      upper = c(pi = 1 - .Machine$double.eps, shape = Inf)
    ))
  }

  if (inherits(model, "logarithmicModel")) {
    start = if ("start" %in% names(dots)) {
      dots$start
    } else if ("pi" %in% names(dots)) {
      dots$pi
    } else {
      0.5
    }
    validateLogarithmicPi(start, "start")
    return(list(
      start = c(pi = start),
      lower = c(pi = sqrt(.Machine$double.eps)),
      upper = c(pi = 1 - sqrt(.Machine$double.eps))
    ))
  }

  stop(
    "weighted MLE fitting currently supports zeta, ziz, and logarithmic models",
    call. = FALSE
  )
}

#' Fit a built-in model using arbitrary positive real observation weights
#'
#' This is the internal fitting boundary required by Rubin Bayesian Bootstrap.
#' It retains the observed `psData` object unchanged, substitutes real-valued
#' weights only in a private survey copy, and evaluates the established model
#' log-likelihood and probability contracts without expanding observations.
#'
#' @param x An object of class `psData` containing the observed survey.
#' @param model A built-in model descriptor returned by [zetaModel()],
#'   [zizModel()], or [logarithmicModel()].
#' @param weights Positive finite fitting weights aligned with occupied support
#'   rows in `x`.
#' @param nterms Number of fitted P/S probability terms to retain.
#' @param ... Optional starting-value controls passed to `weightedMleControl()`.
#' @return An internal `psWeightedMle` object containing fitted parameters,
#'   probabilities, optimiser output, and the supplied weights.
#' @keywords internal
#' @noRd
fitWeightedModel = function(x, model, weights, nterms = 10, ...) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData")
  }
  if (!inherits(model, "psModel")) {
    stop("model must inherit from psModel")
  }
  if (!is.numeric(nterms) || length(nterms) != 1L || !is.finite(nterms) ||
      nterms < 1 || nterms != floor(nterms)) {
    stop("nterms must be one positive integer")
  }

  validateMleObservationSupport(x)
  weights = validateFittingWeights(x, weights)
  weightedData = weightedSurveyData(x, weights)
  control = weightedMleControl(model, ...)
  parameterNames = modelParameterNames(model)

  objective = function(parameterValues) {
    names(parameterValues) = parameterNames
    value = tryCatch(
      modelLogLikelihood(
        model = model,
        parameters = as.list(parameterValues),
        data = weightedData
      ),
      error = function(error) NA_real_
    )
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
      return(.Machine$double.xmax^0.5)
    }
    -value
  }

  fittedModel = optim(
    par = control$start,
    fn = objective,
    method = "L-BFGS-B",
    lower = control$lower,
    upper = control$upper,
    hessian = TRUE
  )
  if (!identical(fittedModel$convergence, 0L)) {
    stop(
      "weighted MLE optimisation failed for model '",
      model$model,
      "': ",
      fittedModel$message,
      call. = FALSE
    )
  }

  parameters = as.numeric(fittedModel$par)
  names(parameters) = parameterNames
  probabilityIndices = posteriorProbabilityIndices(x$type, as.integer(nterms))
  probabilityMatrix = modelProbabilities(
    model = model,
    parameters = as.list(parameters),
    n = probabilityIndices,
    type = x$type
  )
  fittedProbabilities = as.numeric(probabilityMatrix[1L, ])
  names(fittedProbabilities) = colnames(probabilityMatrix)

  result = list(
    psData = x,
    weights = weights,
    parameters = parameters,
    fitted = fittedProbabilities,
    fit = fittedModel,
    model = model$model,
    modelObject = model
  )
  class(result) = "psWeightedMle"
  result
}
