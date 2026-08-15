#' Construct an internal fitPS statistical model descriptor.
#'
#' @param model Non-empty model identifier.
#' @param parameterNames Character vector of natural parameter names.
#' @param supportedEngines Character vector of supported posterior-engine names.
#' @param subclass Concrete S3 subclass placed before `psModel`.
#' @param mleStart Optional named numeric starting values for generic MLE fitting.
#' @param mleLower Optional named numeric lower bounds for generic MLE fitting.
#' @param mleUpper Optional named numeric upper bounds for generic MLE fitting.
#' @return An internal object inheriting from `psModel`.
#' @keywords internal
#' @noRd
newPsModel = function(model,
                       parameterNames,
                       supportedEngines,
                       subclass,
                       mleStart = NULL,
                       mleLower = NULL,
                       mleUpper = NULL) {
  if (!is.character(model) || length(model) != 1L || !nzchar(model)) {
    stop("model must be one non-empty character value")
  }
  if (!is.character(parameterNames) || length(parameterNames) == 0L ||
      any(!nzchar(parameterNames)) || anyDuplicated(parameterNames)) {
    stop("parameterNames must contain unique non-empty character values")
  }
  if (!is.character(supportedEngines) ||
      any(!nzchar(supportedEngines)) || anyDuplicated(supportedEngines)) {
    stop("supportedEngines must contain unique non-empty character values")
  }
  if (!is.character(subclass) || length(subclass) != 1L || !nzchar(subclass)) {
    stop("subclass must be one non-empty character value")
  }

  result = list(
    model = model,
    parameterNames = parameterNames,
    supportedEngines = supportedEngines,
    mleStart = mleStart,
    mleLower = mleLower,
    mleUpper = mleUpper
  )
  class(result) = c(subclass, "psModel")
  result
}


#' Construct a public fitPS model descriptor
#'
#' Creates a model descriptor that can be defined outside fitPS and passed to
#' [fit()]. External model classes supply model-specific mathematics through the
#' exported S3 generics documented here. For simple maximum-likelihood models,
#' fitPS can own the numerical optimisation when starting values and parameter
#' bounds are supplied.
#'
#' @param model Non-empty model identifier stored in fitted objects.
#' @param parameterNames Unique natural-scale parameter names.
#' @param subclass Concrete S3 subclass for external method dispatch.
#' @param supportedEngines Character vector of supported Bayesian posterior
#'   engines. MLE-only models may use `character()`.
#' @param mleStart Optional named numeric starting values for generic MLE.
#' @param mleLower Optional named numeric lower bounds. Missing bounds default to
#'   `-Inf`.
#' @param mleUpper Optional named numeric upper bounds. Missing bounds default to
#'   `Inf`.
#' @return An object inheriting from `psModel` and `subclass`.
#' @export
#'
#' @examples
#' model = psModel(
#'   model = "example",
#'   parameterNames = "theta",
#'   subclass = "exampleModel",
#'   mleStart = c(theta = 1),
#'   mleLower = c(theta = 0)
#' )
psModel = function(model,
                    parameterNames,
                    subclass,
                    supportedEngines = character(),
                    mleStart = NULL,
                    mleLower = NULL,
                    mleUpper = NULL) {
  newPsModel(
    model = model,
    parameterNames = parameterNames,
    supportedEngines = supportedEngines,
    subclass = subclass,
    mleStart = mleStart,
    mleLower = mleLower,
    mleUpper = mleUpper
  )
}

#' Construct built-in fitPS model descriptors
#'
#' These constructors create model objects for the built-in fitPS distributions.
#' They can be supplied directly to [fit()]. The lower-level third-party model
#' constructor and extension generics are introduced separately after the public
#' contract has been validated.
#'
#' @return An object inheriting from `psModel` and the concrete built-in model
#'   class.
#' @export
zetaModel = function() {
  newPsModel(
    model = "zeta",
    parameterNames = "shape",
    supportedEngines = c("numerical", "mcmc"),
    subclass = "zetaModel"
  )
}

#' @rdname zetaModel
#' @export
zizModel = function() {
  newPsModel(
    model = "ziz",
    parameterNames = c("pi", "shape"),
    supportedEngines = c("numerical", "mcmc", "laplace", "importance"),
    subclass = "zizModel"
  )
}

#' @rdname zetaModel
#' @export
logarithmicModel = function() {
  newPsModel(
    model = "logarithmic",
    parameterNames = "pi",
    supportedEngines = c("numerical", "mcmc"),
    subclass = "logarithmicModel"
  )
}

#' Return the natural parameter names declared by a fitPS model.
#'
#' @param model A `psModel` object.
#' @param ... Additional arguments reserved for model methods.
#' @return Character vector of parameter names.
#' @export
modelParameterNames = function(model, ...) {
  UseMethod("modelParameterNames")
}

#' @rdname modelParameterNames
#' @keywords internal
#' @exportS3Method modelParameterNames psModel
#' @noRd
modelParameterNames.psModel = function(model, ...) {
  model$parameterNames
}

#' Return posterior engines supported by a fitPS model.
#'
#' @param model A `psModel` object.
#' @param ... Additional arguments reserved for model methods.
#' @return Character vector of posterior-engine identifiers.
#' @export
supportedPosteriorEngines = function(model, ...) {
  UseMethod("supportedPosteriorEngines")
}

#' @rdname supportedPosteriorEngines
#' @keywords internal
#' @exportS3Method supportedPosteriorEngines psModel
#' @noRd
supportedPosteriorEngines.psModel = function(model, ...) {
  model$supportedEngines
}

#' Test whether a model supports a posterior engine.
#'
#' @param model An internal `psModel` object.
#' @param engine A `psPosteriorEngine` object or one engine name.
#' @return Logical scalar.
#' @keywords internal
#' @noRd
supportsPosteriorEngine = function(model, engine) {
  if (!inherits(model, "psModel")) {
    stop("model must inherit from psModel")
  }

  method = if (inherits(engine, "psPosteriorEngine")) {
    posteriorEngineName(engine)
  } else {
    if (!is.character(engine) || length(engine) != 1L || !nzchar(engine)) {
      stop("engine must be a psPosteriorEngine or one non-empty character value")
    }
    engine
  }

  method %in% supportedPosteriorEngines(model)
}

#' Convert fitPS data to the observation support required by a model.
#'
#' @param model A `psModel` object.
#' @param x An object of class `psData`.
#' @param ... Additional arguments reserved for model methods.
#' @return Model-specific observation data.
#' @export
modelObservationData = function(model, x, ...) {
  UseMethod("modelObservationData")
}

#' @rdname modelObservationData
#' @keywords internal
#' @exportS3Method modelObservationData psModel
#' @noRd
modelObservationData.psModel = function(model, x, ...) {
  psObservationData(x)
}

#' Return generic maximum-likelihood controls for a model
#'
#' External models that use fitPS-owned numerical optimisation may store their
#' starting values and bounds in [psModel()] or override this generic.
#'
#' @param model A `psModel` object.
#' @param x An object of class `psData`.
#' @param ... Additional model-specific fitting controls.
#' @return A list with named numeric `start`, `lower`, and `upper` components.
#' @export
modelMleControl = function(model, x, ...) {
  UseMethod("modelMleControl")
}

#' @rdname modelMleControl
#' @exportS3Method modelMleControl psModel
modelMleControl.psModel = function(model, x, ...) {
  parameterNames = modelParameterNames(model)
  start = model$mleStart
  if (is.null(start)) {
    stop(
      "generic MLE fitting requires starting values; supply mleStart in ",
      "psModel() or define modelMleControl() for class '",
      class(model)[1L],
      "'",
      call. = FALSE
    )
  }

  validateMleVector = function(value, name, allowInfinite = FALSE) {
    if (!is.numeric(value) || length(value) != length(parameterNames) ||
        is.null(names(value)) || !setequal(names(value), parameterNames)) {
      stop(name, " must be a named numeric vector matching parameterNames")
    }
    value = value[parameterNames]
    if (allowInfinite) {
      if (any(is.na(value))) {
        stop(name, " must not contain missing values")
      }
    } else if (any(!is.finite(value))) {
      stop(name, " must contain finite values")
    }
    value
  }

  start = validateMleVector(start, "mleStart")
  lower = if (is.null(model$mleLower)) {
    rep(-Inf, length(parameterNames))
  } else {
    validateMleVector(model$mleLower, "mleLower", allowInfinite = TRUE)
  }
  upper = if (is.null(model$mleUpper)) {
    rep(Inf, length(parameterNames))
  } else {
    validateMleVector(model$mleUpper, "mleUpper", allowInfinite = TRUE)
  }
  names(lower) = parameterNames
  names(upper) = parameterNames

  if (any(lower >= upper)) {
    stop("each MLE lower bound must be less than its upper bound")
  }
  if (any(start <= lower | start >= upper)) {
    stop("each MLE starting value must lie strictly inside its bounds")
  }

  list(start = start, lower = lower, upper = upper)
}

#' Evaluate fitted P/S probabilities through a model descriptor.
#'
#' @param model A `psModel` object.
#' @param parameters Named model parameters.
#' @param n Requested P/S probability indices.
#' @param type Survey type, either `"P"` or `"S"`.
#' @param ... Additional arguments reserved for model methods.
#' @return Numeric matrix of fitted P/S probabilities.
#' @export
modelProbabilities = function(model, parameters, n, type, ...) {
  UseMethod("modelProbabilities")
}

#' Extract one named model parameter from an internal parameter container.
#'
#' @param parameters Named list or named atomic vector.
#' @param name Required parameter name.
#' @return The named parameter value.
#' @keywords internal
#' @noRd
modelParameter = function(parameters, name) {
  if (is.null(names(parameters)) || !name %in% names(parameters)) {
    stop("parameters must contain a named ", name, " component")
  }
  parameters[[name]]
}

#' @rdname modelProbabilities
#' @keywords internal
#' @exportS3Method modelProbabilities zetaModel
#' @noRd
modelProbabilities.zetaModel = function(model, parameters, n, type, ...) {
  zetaProbabilities(
    shape = modelParameter(parameters, "shape"),
    n = n,
    type = type
  )
}

#' @rdname modelProbabilities
#' @keywords internal
#' @exportS3Method modelProbabilities zizModel
#' @noRd
modelProbabilities.zizModel = function(model, parameters, n, type, ...) {
  zizProbabilities(
    pi = modelParameter(parameters, "pi"),
    shape = modelParameter(parameters, "shape"),
    n = n,
    type = type
  )
}

#' Evaluate a model log likelihood through the shared model protocol.
#'
#' Concrete likelihood methods implement the model-specific mathematics behind
#' the shared interface. The base method fails explicitly so an unsupported
#' model cannot silently use the wrong likelihood.
#'
#' @param model A `psModel` object.
#' @param parameters Named model parameters.
#' @param data Model-specific observation data.
#' @param ... Additional model-specific inputs.
#' @return A numeric log likelihood when implemented by a concrete model.
#' @export
modelLogLikelihood = function(model, parameters, data, ...) {
  UseMethod("modelLogLikelihood")
}

#' @rdname modelLogLikelihood
#' @keywords internal
#' @exportS3Method modelLogLikelihood psModel
#' @noRd
modelLogLikelihood.psModel = function(model, parameters, data, ...) {
  stop(
    "modelLogLikelihood is not yet implemented for model '",
    model$model,
    "'",
    call. = FALSE
  )
}

#' Evaluate the plain zeta log likelihood.
#'
#' @param model A `zetaModel` descriptor.
#' @param parameters Named parameters containing scalar `shape`.
#' @param data An object of class `psData`.
#' @param ... Additional arguments reserved for future zeta likelihood controls.
#' @return Scalar log likelihood.
#' @keywords internal
#' @exportS3Method modelLogLikelihood zetaModel
#' @noRd
modelLogLikelihood.zetaModel = function(model, parameters, data, ...) {
  if (!is(data, "psData")) {
    stop("data must be an object of class psData")
  }

  shape = modelParameter(parameters, "shape")
  if (!is.numeric(shape) || length(shape) != 1L || !is.finite(shape)) {
    stop("shape must be one finite numeric value")
  }
  validateZetaShape(shape)

  observations = modelObservationData(model, data)
  sum(data$data$rn * dzetaStandard(observations, shape = shape, log = TRUE))
}

#' Evaluate the zero-inflated zeta log likelihood.
#'
#' @param model A `zizModel` descriptor.
#' @param parameters Named parameters containing scalar `pi` and `shape`.
#' @param data An object of class `psData`.
#' @param ... Additional arguments reserved for future ZIZ likelihood controls.
#' @return Scalar log likelihood.
#' @keywords internal
#' @exportS3Method modelLogLikelihood zizModel
#' @noRd
modelLogLikelihood.zizModel = function(model, parameters, data, ...) {
  if (!is(data, "psData")) {
    stop("data must be an object of class psData")
  }

  pi = modelParameter(parameters, "pi")
  shape = modelParameter(parameters, "shape")

  if (!is.numeric(pi) || length(pi) != 1L || !is.finite(pi)) {
    stop("pi must be one finite numeric value")
  }
  if (!is.numeric(shape) || length(shape) != 1L || !is.finite(shape)) {
    stop("shape must be one finite numeric value")
  }

  observations = modelObservationData(model, data)
  zizLogLikelihood(
    obsData = observations,
    counts = data$data$rn,
    pi = pi,
    shape = shape
  )
}


#' @rdname modelProbabilities
#' @keywords internal
#' @exportS3Method modelProbabilities logarithmicModel
#' @noRd
modelProbabilities.logarithmicModel = function(model, parameters, n, type, ...) {
  logarithmicProbabilities(
    pi = modelParameter(parameters, "pi"),
    n = n,
    type = type
  )
}

#' Evaluate the logarithmic model log likelihood.
#'
#' @param model A `logarithmicModel` descriptor.
#' @param parameters Named parameters containing scalar `pi`.
#' @param data An object of class `psData`.
#' @param ... Additional arguments reserved for future logarithmic controls.
#' @return Scalar log likelihood.
#' @keywords internal
#' @exportS3Method modelLogLikelihood logarithmicModel
#' @noRd
modelLogLikelihood.logarithmicModel = function(model, parameters, data, ...) {
  if (!is(data, "psData")) {
    stop("data must be an object of class psData")
  }

  pi = modelParameter(parameters, "pi")
  validateLogarithmicPi(pi)

  observations = modelObservationData(model, data)
  sum(data$data$rn * dlog(observations, shape = pi, log = TRUE))
}
