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
#' External methods are responsible for mapping the P/S survey labels onto the
#' model's natural support. For example, a zero-based distribution may use P0,
#' P1, ... directly but map S1, S2, ... to support values 0, 1, ... so the
#' probability sequence is shifted without truncation or renormalisation.
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
#' @param n Requested P/S probability indices. These are survey labels, not
#'   necessarily the model's natural support values; concrete methods should
#'   perform any required support shift before evaluating probabilities.
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
