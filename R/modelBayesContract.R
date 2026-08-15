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
#'   engines may require a named numeric `start` component. One-dimensional
#'   numerical fitting also uses named numeric `lower` and `upper` components.
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

#' Transform natural model parameters to unconstrained coordinates
#'
#' Generic MCMC proposes parameter values on an unconstrained Euclidean scale.
#' This avoids proposals that repeatedly leave a model's natural parameter
#' space, such as negative values for a positive rate or values outside
#' `(0, 1)` for a probability. External models with constrained parameters
#' implement this generic to map their named natural-scale parameter vector to
#' a named vector whose components may vary over the whole real line.
#'
#' The transformation must be one-to-one with [modelFromUnconstrained()] over
#' the parameter region used for Bayesian fitting, preserve the names and order
#' returned by [modelParameterNames()], and return finite values for valid
#' interior natural-scale parameters. The default `psModel` method is the
#' identity transformation, which is appropriate when every natural parameter
#' is already unconstrained.
#'
#' Common choices include `log(theta)` for a positive parameter `theta` and
#' `qlogis(p)` for a probability `p`. If several parameters are transformed,
#' the returned vector contains one unconstrained coordinate for each natural
#' parameter.
#'
#' @param model A `psModel` object.
#' @param parameters Named numeric model parameters on their natural scale.
#' @param ... Additional model-specific transformation inputs.
#' @return A named numeric vector of unconstrained coordinates, with names and
#'   order matching `modelParameterNames(model)`.
#' @seealso [modelFromUnconstrained()], [modelLogJacobian()]
#' @export
modelToUnconstrained = function(model, parameters, ...) {
  UseMethod("modelToUnconstrained")
}

#' @rdname modelToUnconstrained
#' @exportS3Method modelToUnconstrained psModel
modelToUnconstrained.psModel = function(model, parameters, ...) {
  validateModelParameterVector(model, parameters, "parameters")
}

#' Transform unconstrained coordinates to natural model parameters
#'
#' This is the inverse of [modelToUnconstrained()]. Generic MCMC calls it for
#' every proposed unconstrained state before evaluating the likelihood and
#' prior, both of which remain defined on the model's natural parameter scale.
#'
#' The transformation must be one-to-one with [modelToUnconstrained()] over the
#' Bayesian parameter region, preserve the names and order returned by
#' [modelParameterNames()], and return valid natural-scale parameters whenever
#' its input is finite and valid for the chosen transformation. The default
#' `psModel` method is the identity transformation.
#'
#' For example, if `modelToUnconstrained()` uses `log(theta)` for a positive
#' parameter, this method uses `exp(z)`. If it uses `qlogis(p)` for a
#' probability, this method uses `plogis(z)`.
#'
#' @param model A `psModel` object.
#' @param unconstrained Named numeric coordinates on the unconstrained scale.
#' @param ... Additional model-specific transformation inputs.
#' @return A named numeric vector of natural-scale model parameters, with names
#'   and order matching `modelParameterNames(model)`.
#' @seealso [modelToUnconstrained()], [modelLogJacobian()]
#' @export
modelFromUnconstrained = function(model, unconstrained, ...) {
  UseMethod("modelFromUnconstrained")
}

#' @rdname modelFromUnconstrained
#' @exportS3Method modelFromUnconstrained psModel
modelFromUnconstrained.psModel = function(model, unconstrained, ...) {
  validateModelParameterVector(model, unconstrained, "unconstrained")
}

#' Evaluate the log-Jacobian for the unconstrained parameter transformation
#'
#' When MCMC samples an unconstrained vector `z` but the model prior is defined
#' for natural parameters `theta = modelFromUnconstrained(model, z)`, the
#' posterior density on the sampled scale contains the change-of-variables
#' factor
#'
#' `log |det(d theta / d z)|`.
#'
#' This generic returns exactly that quantity: the log absolute determinant of
#' the Jacobian of the transformation **from unconstrained coordinates back to
#' natural parameters**. It is not the Jacobian of
#' [modelToUnconstrained()]. Generic MCMC adds this value to the natural-scale
#' log likelihood and log prior.
#'
#' For an identity transformation the value is zero. For `theta = exp(z)` it is
#' `z`. For `p = plogis(z)` it is `log(p) + log1p(-p)`. For several independent
#' component-wise transformations, return the sum of their log-Jacobian terms;
#' for a genuinely coupled transformation, return the log absolute determinant
#' of the full inverse-transformation Jacobian matrix.
#'
#' @param model A `psModel` object.
#' @param unconstrained Named numeric coordinates on the unconstrained scale.
#' @param ... Additional model-specific transformation inputs.
#' @return One finite numeric value equal to
#'   `log |det(d theta / d z)|` for the inverse transformation.
#' @seealso [modelToUnconstrained()], [modelFromUnconstrained()]
#' @export
modelLogJacobian = function(model, unconstrained, ...) {
  UseMethod("modelLogJacobian")
}

#' @rdname modelLogJacobian
#' @exportS3Method modelLogJacobian psModel
modelLogJacobian.psModel = function(model, unconstrained, ...) {
  validateModelParameterVector(model, unconstrained, "unconstrained")
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

  list(
    start = c(shape = start),
    lower = c(shape = prior$range[1L]),
    upper = c(shape = prior$range[2L])
  )
}

#' @rdname modelToUnconstrained
#' @exportS3Method modelToUnconstrained zetaModel
modelToUnconstrained.zetaModel = function(model, parameters, ...) {
  parameters = validateModelParameterVector(model, parameters, "parameters")
  c(shape = shapeToTau(parameters[["shape"]]))
}

#' @rdname modelFromUnconstrained
#' @exportS3Method modelFromUnconstrained zetaModel
modelFromUnconstrained.zetaModel = function(model, unconstrained, ...) {
  unconstrained = validateModelParameterVector(model, unconstrained, "unconstrained")
  c(shape = tauToShape(unconstrained[["shape"]]))
}

#' @rdname modelLogJacobian
#' @exportS3Method modelLogJacobian zetaModel
modelLogJacobian.zetaModel = function(model, unconstrained, ...) {
  unconstrained = validateModelParameterVector(model, unconstrained, "unconstrained")
  unname(unconstrained[["shape"]])
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

  list(
    start = c(pi = start),
    lower = c(pi = prior$range[1L]),
    upper = c(pi = prior$range[2L])
  )
}

#' @rdname modelToUnconstrained
#' @exportS3Method modelToUnconstrained logarithmicModel
modelToUnconstrained.logarithmicModel = function(model, parameters, ...) {
  parameters = validateModelParameterVector(model, parameters, "parameters")
  c(pi = logitPi(parameters[["pi"]]))
}

#' @rdname modelFromUnconstrained
#' @exportS3Method modelFromUnconstrained logarithmicModel
modelFromUnconstrained.logarithmicModel = function(model, unconstrained, ...) {
  unconstrained = validateModelParameterVector(model, unconstrained, "unconstrained")
  c(pi = invLogitPi(unconstrained[["pi"]]))
}

#' @rdname modelLogJacobian
#' @exportS3Method modelLogJacobian logarithmicModel
modelLogJacobian.logarithmicModel = function(model, unconstrained, ...) {
  unconstrained = validateModelParameterVector(model, unconstrained, "unconstrained")
  pi = invLogitPi(unconstrained[["pi"]])
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

  boundaryInset = sqrt(.Machine$double.eps)
  list(
    start = start,
    lower = c(pi = boundaryInset, shape = prior$range[[1L]]),
    upper = c(pi = 1 - boundaryInset, shape = prior$range[[2L]])
  )
}

#' @rdname modelToUnconstrained
#' @exportS3Method modelToUnconstrained zizModel
modelToUnconstrained.zizModel = function(model, parameters, ...) {
  parameters = validateModelParameterVector(model, parameters, "parameters")
  unconstrained = zizThetaToWorking(parameters)
  names(unconstrained) = modelParameterNames(model)
  unconstrained
}

#' @rdname modelFromUnconstrained
#' @exportS3Method modelFromUnconstrained zizModel
modelFromUnconstrained.zizModel = function(model, unconstrained, ...) {
  unconstrained = validateModelParameterVector(model, unconstrained, "unconstrained")
  parameters = zizWorkingToTheta(unconstrained)
  parameters[modelParameterNames(model)]
}

#' @rdname modelLogJacobian
#' @exportS3Method modelLogJacobian zizModel
modelLogJacobian.zizModel = function(model, unconstrained, ...) {
  unconstrained = validateModelParameterVector(model, unconstrained, "unconstrained")
  zizWorkingLogJacobian(unconstrained)
}
