#' Construct an internal fitPS statistical model descriptor.
#'
#' @param model Non-empty model identifier.
#' @param parameterNames Character vector of natural parameter names.
#' @param supportedEngines Character vector of supported posterior-engine names.
#' @param subclass Concrete S3 subclass placed before `psModel`.
#' @return An internal object inheriting from `psModel`.
#' @keywords internal
#' @noRd
newPsModel = function(model,
                       parameterNames,
                       supportedEngines,
                       subclass) {
  if (!is.character(model) || length(model) != 1L || !nzchar(model)) {
    stop("model must be one non-empty character value")
  }
  if (!is.character(parameterNames) || length(parameterNames) == 0L ||
      any(!nzchar(parameterNames)) || anyDuplicated(parameterNames)) {
    stop("parameterNames must contain unique non-empty character values")
  }
  if (!is.character(supportedEngines) || length(supportedEngines) == 0L ||
      any(!nzchar(supportedEngines)) || anyDuplicated(supportedEngines)) {
    stop("supportedEngines must contain unique non-empty character values")
  }
  if (!is.character(subclass) || length(subclass) != 1L || !nzchar(subclass)) {
    stop("subclass must be one non-empty character value")
  }

  result = list(
    model = model,
    parameterNames = parameterNames,
    supportedEngines = supportedEngines
  )
  class(result) = c(subclass, "psModel")
  result
}

#' Construct the plain zeta model descriptor.
#'
#' @return An internal object inheriting from `zetaModel` and `psModel`.
#' @keywords internal
#' @noRd
zetaModel = function() {
  newPsModel(
    model = "zeta",
    parameterNames = "shape",
    supportedEngines = c("numerical", "mcmc"),
    subclass = "zetaModel"
  )
}

#' Construct the zero-inflated zeta model descriptor.
#'
#' @return An internal object inheriting from `zizModel` and `psModel`.
#' @keywords internal
#' @noRd
zizModel = function() {
  newPsModel(
    model = "ziz",
    parameterNames = c("pi", "shape"),
    supportedEngines = c("numerical", "mcmc", "laplace", "importance"),
    subclass = "zizModel"
  )
}

#' Return the natural parameter names declared by a fitPS model.
#'
#' @param model An internal `psModel` object.
#' @param ... Additional arguments reserved for model methods.
#' @return Character vector of parameter names.
#' @keywords internal
#' @noRd
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
#' @param model An internal `psModel` object.
#' @param ... Additional arguments reserved for model methods.
#' @return Character vector of posterior-engine identifiers.
#' @keywords internal
#' @noRd
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
#' @param model An internal `psModel` object.
#' @param x An object of class `psData`.
#' @param ... Additional arguments reserved for model methods.
#' @return Model-specific observation data.
#' @keywords internal
#' @noRd
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

#' Evaluate fitted P/S probabilities through a model descriptor.
#'
#' @param model An internal `psModel` object.
#' @param parameters Named model parameters.
#' @param n Requested P/S probability indices.
#' @param type Survey type, either `"P"` or `"S"`.
#' @param ... Additional arguments reserved for model methods.
#' @return Numeric matrix of fitted P/S probabilities.
#' @keywords internal
#' @noRd
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

#' Evaluate a model log likelihood through the Stage 6 model protocol.
#'
#' Concrete likelihood methods are introduced when each fitting path migrates
#' to the new architecture. The base method fails explicitly so an unsupported
#' or not-yet-migrated path cannot silently use the wrong mathematics.
#'
#' @param model An internal `psModel` object.
#' @param parameters Named model parameters.
#' @param data Model-specific observation data.
#' @param ... Additional model-specific inputs.
#' @return A numeric log likelihood when implemented by a concrete model.
#' @keywords internal
#' @noRd
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
