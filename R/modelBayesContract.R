#' Evaluate a model prior on the natural parameter scale
#'
#' External Bayesian models implement this generic to provide the complete
#' log-prior density for a named natural-scale parameter vector. The prior
#' object is deliberately model-defined so multi-parameter and correlated
#' priors do not have to use the legacy one-dimensional `psPrior` structure.
#'
#' @param model A `psModel` object.
#' @param parameters Named natural-scale model parameters.
#' @param prior A model-specific prior specification.
#' @param shape1,shape2 Positive beta-prior shape parameters used by models that
#'   place a beta prior on a probability parameter. Ignored by models that do
#'   not use these arguments.
#' @param ... Additional model-specific prior inputs.
#' @return One numeric log-prior density value, or `-Inf` outside prior support.
#' @export
modelLogPrior = function(model, parameters, prior, ...) {
  UseMethod("modelLogPrior")
}

#' @rdname modelLogPrior
#' @exportS3Method modelLogPrior psModel
modelLogPrior.psModel = function(model, parameters, prior, ...) {
  stop(
    "modelLogPrior is not yet implemented for model '",
    model$model,
    "'",
    call. = FALSE
  )
}

#' Return Bayesian fitting controls supplied by a model
#'
#' This generic exposes model-owned Bayesian controls such as natural-scale
#' starting values. Engine-owned controls such as iteration counts, tolerances,
#' and random seeds remain arguments to the posterior engine.
#'
#' @param model A `psModel` object.
#' @param x An object of class `psData`.
#' @param engine A posterior-engine object or engine name.
#' @param prior A model-specific prior specification.
#' @param start Optional named numeric natural-scale starting values used by
#'   models that support caller-supplied Bayesian starts. Ignored by models that
#'   do not use this argument.
#' @param ... Additional model-specific control inputs.
#' @return A named list containing model-specific Bayesian controls. Generic
#'   engines may require a named numeric `start` component.
#' @export
modelBayesControl = function(model, x, engine, prior, ...) {
  UseMethod("modelBayesControl")
}

#' @rdname modelBayesControl
#' @exportS3Method modelBayesControl psModel
modelBayesControl.psModel = function(model, x, engine, prior, ...) {
  stop(
    "generic Bayesian fitting requires model controls; define ",
    "modelBayesControl() for class '",
    class(model)[1L],
    "'",
    call. = FALSE
  )
}

validateModelParameterVector = function(model, value, name) {
  parameterNames = modelParameterNames(model)

  if (!is.numeric(value) || length(value) != length(parameterNames) ||
      is.null(names(value)) || !setequal(names(value), parameterNames)) {
    stop(name, " must be a named numeric vector matching modelParameterNames(model)")
  }

  value = value[parameterNames]
  if (any(!is.finite(value))) {
    stop(name, " must contain finite values")
  }

  value
}

#' Transform natural model parameters to working coordinates
#'
#' The default `psModel` method is the identity transformation. Models with
#' constrained parameters may override it to expose unconstrained coordinates
#' suitable for generic Bayesian engines.
#'
#' @param model A `psModel` object.
#' @param parameters Named natural-scale model parameters.
#' @param ... Additional model-specific transformation inputs.
#' @return A named numeric working-scale parameter vector.
#' @export
modelToWorking = function(model, parameters, ...) {
  UseMethod("modelToWorking")
}

#' @rdname modelToWorking
#' @exportS3Method modelToWorking psModel
modelToWorking.psModel = function(model, parameters, ...) {
  validateModelParameterVector(model, parameters, "parameters")
}

#' Transform working coordinates to natural model parameters
#'
#' The default `psModel` method is the identity transformation.
#'
#' @param model A `psModel` object.
#' @param working Named working-scale model parameters.
#' @param ... Additional model-specific transformation inputs.
#' @return A named numeric natural-scale parameter vector.
#' @export
modelFromWorking = function(model, working, ...) {
  UseMethod("modelFromWorking")
}

#' @rdname modelFromWorking
#' @exportS3Method modelFromWorking psModel
modelFromWorking.psModel = function(model, working, ...) {
  validateModelParameterVector(model, working, "working")
}

#' Evaluate the inverse-transformation log-Jacobian
#'
#' Returns the log absolute determinant of the Jacobian for the transformation
#' from working coordinates to natural model parameters. The default identity
#' transformation therefore returns zero.
#'
#' @param model A `psModel` object.
#' @param working Named working-scale model parameters.
#' @param ... Additional model-specific transformation inputs.
#' @return One finite numeric log absolute Jacobian value.
#' @export
modelWorkingLogJacobian = function(model, working, ...) {
  UseMethod("modelWorkingLogJacobian")
}

#' @rdname modelWorkingLogJacobian
#' @exportS3Method modelWorkingLogJacobian psModel
modelWorkingLogJacobian.psModel = function(model, working, ...) {
  validateModelParameterVector(model, working, "working")
  0
}

#' @rdname modelLogPrior
#' @exportS3Method modelLogPrior zetaModel
modelLogPrior.zetaModel = function(model, parameters, prior, ...) {
  shape = modelParameter(parameters, "shape")
  validateBayesPrior(prior)
  validateZetaPriorRange(prior$range)
  prior$logd(shape)
}

#' @rdname modelBayesControl
#' @exportS3Method modelBayesControl zetaModel
modelBayesControl.zetaModel = function(model, x, engine, prior, ...) {
  validateBayesPrior(prior)
  validateZetaPriorRange(prior$range)

  start = 2
  if (!inRange(start, prior$range)) {
    start = mean(prior$range)
  }

  list(start = c(shape = start))
}

#' @rdname modelToWorking
#' @exportS3Method modelToWorking zetaModel
modelToWorking.zetaModel = function(model, parameters, ...) {
  parameters = validateModelParameterVector(model, parameters, "parameters")
  c(shape = shapeToTau(parameters[["shape"]]))
}

#' @rdname modelFromWorking
#' @exportS3Method modelFromWorking zetaModel
modelFromWorking.zetaModel = function(model, working, ...) {
  working = validateModelParameterVector(model, working, "working")
  c(shape = tauToShape(working[["shape"]]))
}

#' @rdname modelWorkingLogJacobian
#' @exportS3Method modelWorkingLogJacobian zetaModel
modelWorkingLogJacobian.zetaModel = function(model, working, ...) {
  working = validateModelParameterVector(model, working, "working")
  unname(working[["shape"]])
}

#' @rdname modelLogPrior
#' @exportS3Method modelLogPrior logarithmicModel
modelLogPrior.logarithmicModel = function(model, parameters, prior, ...) {
  pi = modelParameter(parameters, "pi")
  validateBayesPrior(prior)
  validateLogarithmicPriorRange(prior$range)
  prior$logd(pi)
}

#' @rdname modelBayesControl
#' @exportS3Method modelBayesControl logarithmicModel
modelBayesControl.logarithmicModel = function(model, x, engine, prior, ...) {
  validateBayesPrior(prior)
  validateLogarithmicPriorRange(prior$range)

  start = 0.5
  if (!inRange(start, prior$range)) {
    start = mean(prior$range)
  }

  list(start = c(pi = start))
}

#' @rdname modelToWorking
#' @exportS3Method modelToWorking logarithmicModel
modelToWorking.logarithmicModel = function(model, parameters, ...) {
  parameters = validateModelParameterVector(model, parameters, "parameters")
  c(pi = logitPi(parameters[["pi"]]))
}

#' @rdname modelFromWorking
#' @exportS3Method modelFromWorking logarithmicModel
modelFromWorking.logarithmicModel = function(model, working, ...) {
  working = validateModelParameterVector(model, working, "working")
  c(pi = invLogitPi(working[["pi"]]))
}

#' @rdname modelWorkingLogJacobian
#' @exportS3Method modelWorkingLogJacobian logarithmicModel
modelWorkingLogJacobian.logarithmicModel = function(model, working, ...) {
  working = validateModelParameterVector(model, working, "working")
  pi = invLogitPi(working[["pi"]])
  unname(log(pi) + log1p(-pi))
}

#' @rdname modelLogPrior
#' @exportS3Method modelLogPrior zizModel
modelLogPrior.zizModel = function(model,
                                  parameters,
                                  prior,
                                  shape1 = 1,
                                  shape2 = 1,
                                  ...) {
  pi = modelParameter(parameters, "pi")
  shape = modelParameter(parameters, "shape")

  validateBayesPrior(prior)
  validateZetaPriorRange(prior$range)
  if (!is.numeric(shape1) || length(shape1) != 1L || !is.finite(shape1) || shape1 <= 0) {
    stop("shape1 must be one positive finite numeric value")
  }
  if (!is.numeric(shape2) || length(shape2) != 1L || !is.finite(shape2) || shape2 <= 0) {
    stop("shape2 must be one positive finite numeric value")
  }

  dbeta(pi, shape1 = shape1, shape2 = shape2, log = TRUE) + prior$logd(shape)
}

#' @rdname modelBayesControl
#' @exportS3Method modelBayesControl zizModel
modelBayesControl.zizModel = function(model,
                                      x,
                                      engine,
                                      prior,
                                      start = c(pi = 0.5, shape = 2),
                                      ...) {
  validateBayesPrior(prior)
  validateZetaPriorRange(prior$range)
  start = validateModelParameterVector(model, start, "start")

  if (start[["pi"]] <= 0 || start[["pi"]] >= 1) {
    stop("start pi must lie strictly between 0 and 1")
  }
  if (!inRange(start[["shape"]], prior$range)) {
    start[["shape"]] = mean(prior$range)
  }

  list(start = start)
}

#' @rdname modelToWorking
#' @exportS3Method modelToWorking zizModel
modelToWorking.zizModel = function(model, parameters, ...) {
  parameters = validateModelParameterVector(model, parameters, "parameters")
  working = zizThetaToWorking(parameters)
  names(working) = modelParameterNames(model)
  working
}

#' @rdname modelFromWorking
#' @exportS3Method modelFromWorking zizModel
modelFromWorking.zizModel = function(model, working, ...) {
  working = validateModelParameterVector(model, working, "working")
  parameters = zizWorkingToTheta(working)
  parameters[modelParameterNames(model)]
}

#' @rdname modelWorkingLogJacobian
#' @exportS3Method modelWorkingLogJacobian zizModel
modelWorkingLogJacobian.zizModel = function(model, working, ...) {
  working = validateModelParameterVector(model, working, "working")
  zizWorkingLogJacobian(working)
}
