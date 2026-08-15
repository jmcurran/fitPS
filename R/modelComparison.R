#' Resolve the model descriptor for a fitted fitPS object.
#'
#' @param object An object of class `psFit`.
#' @return An internal `psModel` descriptor.
#' @keywords internal
#' @noRd
modelFromFit = function(object) {
  if (!is(object, "psFit")) {
    stop("object must be an object of class psFit")
  }

  switch(
    object$model,
    zeta = zetaModel(),
    ziz = zizModel(),
    logarithmic = logarithmicModel(),
    log = logarithmicModel(),
    stop("Unknown fitPS model '", object$model, "'", call. = FALSE)
  )
}

#' Extract natural-scale fitted model parameters.
#'
#' @param object An object of class `psFit`.
#' @param model Internal model descriptor for `object`.
#' @return Named numeric vector of model parameters.
#' @keywords internal
#' @noRd
fitModelParameters = function(object, model = modelFromFit(object)) {
  parameterNames = modelParameterNames(model)
  values = vapply(parameterNames, function(parameterName) {
    value = object[[parameterName]]
    if (is.null(value) && !is.null(object$fit$par)) {
      value = object$fit$par[[parameterName]]
    }
    if (is.null(value)) {
      stop("fit does not contain parameter '", parameterName, "'", call. = FALSE)
    }
    as.numeric(value)[1L]
  }, numeric(1L))
  names(values) = parameterNames
  values
}

#' Count free parameters declared by a fitPS model.
#'
#' @param model Internal `psModel` descriptor.
#' @return Integer number of free model parameters.
#' @keywords internal
#' @noRd
modelParameterCount = function(model) {
  length(modelParameterNames(model))
}

#' Evaluate model deviance.
#'
#' The fitPS deviance is `-2 * log L(theta)`. Additive constants independent of
#' the model parameters are omitted, as is standard for likelihood comparison.
#'
#' @param model Internal `psModel` descriptor.
#' @param parameters Named natural-scale model parameters.
#' @param data An object of class `psData`.
#' @return Scalar deviance.
#' @keywords internal
#' @noRd
modelDeviance = function(model, parameters, data) {
  -2 * modelLogLikelihood(model, parameters = parameters, data = data)
}

#' Extract deviance from a fitted fitPS model
#'
#' @param object An object of class `psFit`.
#' @param ... Additional arguments retained for S3 compatibility.
#' @return Scalar model deviance at the fitted representative parameter value.
#' @export
#' @exportS3Method stats::deviance psFit
deviance.psFit = function(object, ...) {
  model = modelFromFit(object)
  modelDeviance(model, fitModelParameters(object, model), object$psData)
}

#' Compute the deviance information criterion for a Bayesian fitPS model
#'
#' DIC is computed as `2 * E[D(theta) | y] - D(E[theta | y])`, where
#' `D(theta) = -2 log L(theta)`. Numerical integration, MCMC draws, and
#' importance weights are taken from the common posterior representation.
#'
#' Laplace fits currently do not retain enough posterior integration information
#' for a stable DIC calculation and therefore fail explicitly rather than using
#' an unrelated approximation silently.
#'
#' @param object A Bayesian object of class `psFit`.
#' @param ... Additional arguments reserved for future controls.
#' @return A scalar numeric DIC value with `pD`, `Dbar`, and `Dhat` attributes.
#' @export
DIC = function(object, ...) {
  if (!is(object, "psFit") || !identical(object$method, "bayes") ||
      is.null(object$posterior)) {
    stop("DIC requires a Bayesian psFit object", call. = FALSE)
  }

  model = modelFromFit(object)
  representation = object$posterior$representation
  dbar = posteriorExpectedDeviance(
    representation = representation,
    model = model,
    data = object$psData
  )
  posteriorMean = object$posterior$parameters$estimate
  names(posteriorMean) = object$posterior$parameters$parameter
  dhat = modelDeviance(model, posteriorMean, object$psData)
  pD = dbar - dhat
  value = dbar + pD
  attr(value, "pD") = pD
  attr(value, "Dbar") = dbar
  attr(value, "Dhat") = dhat
  value
}

#' Compute posterior expected deviance from a posterior representation.
#'
#' @param representation Internal `psPosteriorRepresentation` object.
#' @param model Internal `psModel` descriptor.
#' @param data An object of class `psData`.
#' @param ... Additional representation-specific controls.
#' @return Scalar posterior expected deviance.
#' @keywords internal
#' @noRd
posteriorExpectedDeviance = function(representation, model, data, ...) {
  UseMethod("posteriorExpectedDeviance")
}

#' @rdname posteriorExpectedDeviance
#' @keywords internal
#' @exportS3Method posteriorExpectedDeviance numericalPosteriorRepresentation
#' @noRd
posteriorExpectedDeviance.numericalPosteriorRepresentation = function(representation,
                                                                       model,
                                                                       data,
                                                                       ...) {
  parameterNames = modelParameterNames(model)
  if (length(parameterNames) == 1L &&
      is.function(representation$value$density)) {
    parameterName = parameterNames[1L]
    densityFunction = representation$value$density
    bounds = representation$metadata$bounds
    integral = integrate(function(parameterValue) {
      vapply(parameterValue, function(value) {
        parameters = value
        names(parameters) = parameterName
        modelDeviance(model, parameters, data) * densityFunction(value)
      }, numeric(1L))
    }, lower = bounds[["lower"]], upper = bounds[["upper"]])
    return(integral$value)
  }

  if (inherits(model, "zizModel")) {
    grid = representation$value$posteriorGrid
    deviances = outer(grid$pi, grid$shape, Vectorize(function(pi, shape) {
      modelDeviance(model, c(pi = pi, shape = shape), data)
    }))
    return(sum(deviances * grid$density) * grid$dPi * grid$dShape)
  }

  stop("Numerical DIC is not implemented for model '", model$model, "'", call. = FALSE)
}

#' @rdname posteriorExpectedDeviance
#' @keywords internal
#' @exportS3Method posteriorExpectedDeviance mcmcPosteriorRepresentation
#' @noRd
posteriorExpectedDeviance.mcmcPosteriorRepresentation = function(representation,
                                                                  model,
                                                                  data,
                                                                  ...) {
  chain = representation$value$chain
  if (is.data.frame(chain) || is.matrix(chain)) {
    values = apply(chain, 1L, function(row) {
      row = as.numeric(row)
      names(row) = modelParameterNames(model)
      modelDeviance(model, row, data)
    })
  } else {
    parameterNames = modelParameterNames(model)
    if (length(parameterNames) != 1L) {
      stop("scalar MCMC chains require a one-parameter model")
    }
    values = vapply(chain, function(value) {
      parameters = value
      names(parameters) = parameterNames
      modelDeviance(model, parameters, data)
    }, numeric(1L))
  }
  mean(values)
}

#' @rdname posteriorExpectedDeviance
#' @keywords internal
#' @exportS3Method posteriorExpectedDeviance importancePosteriorRepresentation
#' @noRd
posteriorExpectedDeviance.importancePosteriorRepresentation = function(representation,
                                                                        model,
                                                                        data,
                                                                        ...) {
  samples = representation$value$approximation$samples
  parameterNames = modelParameterNames(model)
  values = apply(samples[parameterNames], 1L, function(row) {
    row = as.numeric(row)
    names(row) = parameterNames
    modelDeviance(model, row, data)
  })
  sum(values * samples$weight) / sum(samples$weight)
}

#' @rdname posteriorExpectedDeviance
#' @keywords internal
#' @exportS3Method posteriorExpectedDeviance laplacePosteriorRepresentation
#' @noRd
posteriorExpectedDeviance.laplacePosteriorRepresentation = function(representation,
                                                                     model,
                                                                     data,
                                                                     ...) {
  stop(
    "DIC is not currently available for Laplace posterior fits because the ",
    "stored Laplace representation does not provide a posterior deviance expectation",
    call. = FALSE
  )
}
