#' Construct an internal posterior-engine descriptor.
#'
#' @param method Non-empty posterior approximation identifier.
#' @param subclass Concrete S3 subclass placed before `psPosteriorEngine`.
#' @return An internal object inheriting from `psPosteriorEngine`.
#' @keywords internal
#' @noRd
newPsPosteriorEngine = function(method, subclass) {
  if (!is.character(method) || length(method) != 1L || !nzchar(method)) {
    stop("method must be one non-empty character value")
  }
  if (!is.character(subclass) || length(subclass) != 1L || !nzchar(subclass)) {
    stop("subclass must be one non-empty character value")
  }

  result = list(method = method)
  class(result) = c(subclass, "psPosteriorEngine")
  result
}

#' Construct the numerical posterior engine descriptor.
#' @return An internal numerical posterior engine.
#' @keywords internal
#' @noRd
numericalPosteriorEngine = function() {
  newPsPosteriorEngine("numerical", "numericalPosteriorEngine")
}

#' Construct the MCMC posterior engine descriptor.
#' @return An internal MCMC posterior engine.
#' @keywords internal
#' @noRd
mcmcPosteriorEngine = function() {
  newPsPosteriorEngine("mcmc", "mcmcPosteriorEngine")
}

#' Construct the Laplace posterior engine descriptor.
#' @return An internal Laplace posterior engine.
#' @keywords internal
#' @noRd
laplacePosteriorEngine = function() {
  newPsPosteriorEngine("laplace", "laplacePosteriorEngine")
}

#' Construct the importance-sampling posterior engine descriptor.
#' @return An internal importance-sampling posterior engine.
#' @keywords internal
#' @noRd
importancePosteriorEngine = function() {
  newPsPosteriorEngine("importance", "importancePosteriorEngine")
}

#' Return the identifier for a posterior engine.
#'
#' @param engine An internal `psPosteriorEngine` object.
#' @param ... Additional arguments reserved for engine methods.
#' @return One posterior approximation identifier.
#' @keywords internal
#' @noRd
posteriorEngineName = function(engine, ...) {
  UseMethod("posteriorEngineName")
}

#' @rdname posteriorEngineName
#' @keywords internal
#' @noRd
posteriorEngineName.psPosteriorEngine = function(engine, ...) {
  engine$method
}

#' Construct an engine-specific posterior representation wrapper.
#'
#' The wrapper gives all posterior representations a small common contract
#' without forcing grids, chains, Gaussian approximations, and weighted
#' samples into one physical data structure.
#'
#' @param engine A `psPosteriorEngine` object.
#' @param value Engine-specific representation payload.
#' @param metadata Optional named list of representation metadata.
#' @return An object inheriting from `psPosteriorRepresentation`.
#' @keywords internal
#' @noRd
newPsPosteriorRepresentation = function(engine, value, metadata = list()) {
  if (!inherits(engine, "psPosteriorEngine")) {
    stop("engine must inherit from psPosteriorEngine")
  }
  if (!is.list(metadata)) {
    stop("metadata must be a list")
  }

  method = posteriorEngineName(engine)
  subclass = paste0(method, "PosteriorRepresentation")
  result = list(
    method = method,
    value = value,
    metadata = metadata
  )
  class(result) = c(subclass, "psPosteriorRepresentation")
  result
}

#' Validate a posterior representation against an optional engine.
#'
#' @param representation An internal `psPosteriorRepresentation` object.
#' @param engine Optional `psPosteriorEngine` expected to own the representation.
#' @return The representation, invisibly.
#' @keywords internal
#' @noRd
validatePosteriorRepresentation = function(representation, engine = NULL) {
  if (!inherits(representation, "psPosteriorRepresentation")) {
    stop("representation must inherit from psPosteriorRepresentation")
  }

  if (!is.null(engine)) {
    if (!inherits(engine, "psPosteriorEngine")) {
      stop("engine must inherit from psPosteriorEngine")
    }
    if (!identical(representation$method, posteriorEngineName(engine))) {
      stop("representation method does not match the posterior engine")
    }
  }

  invisible(representation)
}

#' Fit a posterior through a posterior-engine strategy.
#'
#' @param engine An internal `psPosteriorEngine` object.
#' @param model An internal `psModel` object.
#' @param x An object of class `psData`.
#' @param prior A prior object created by `makePrior()`.
#' @param ... Engine- and model-specific fitting controls.
#' @return An engine-specific posterior representation when implemented.
#' @keywords internal
#' @noRd
fitPosterior = function(engine, model, x, prior, ...) {
  UseMethod("fitPosterior")
}

#' @rdname fitPosterior
#' @keywords internal
#' @noRd
fitPosterior.psPosteriorEngine = function(engine, model, x, prior, ...) {
  validateEngineModelPair(engine, model)
  stop(
    "fitPosterior is not yet implemented for engine '",
    posteriorEngineName(engine),
    "' and model '",
    model$model,
    "'",
    call. = FALSE
  )
}

#' Summarise an engine-specific posterior representation.
#'
#' @param engine An internal `psPosteriorEngine` object.
#' @param model An internal `psModel` object.
#' @param representation A `psPosteriorRepresentation` object.
#' @param ... Engine- and model-specific summary controls.
#' @return A common posterior summary when implemented.
#' @keywords internal
#' @noRd
summarisePosterior = function(engine, model, representation, ...) {
  UseMethod("summarisePosterior")
}

#' @rdname summarisePosterior
#' @keywords internal
#' @noRd
summarisePosterior.psPosteriorEngine = function(engine,
                                                 model,
                                                 representation,
                                                 ...) {
  validateEngineModelPair(engine, model)
  validatePosteriorRepresentation(representation, engine)
  stop(
    "summarisePosterior is not yet implemented for engine '",
    posteriorEngineName(engine),
    "'",
    call. = FALSE
  )
}

#' Extract posterior diagnostics through the engine protocol.
#'
#' @param engine An internal `psPosteriorEngine` object.
#' @param representation A `psPosteriorRepresentation` object.
#' @param ... Engine-specific diagnostic controls.
#' @return Engine diagnostics when implemented.
#' @keywords internal
#' @noRd
posteriorDiagnostics = function(engine, representation, ...) {
  UseMethod("posteriorDiagnostics")
}

#' @rdname posteriorDiagnostics
#' @keywords internal
#' @noRd
posteriorDiagnostics.psPosteriorEngine = function(engine, representation, ...) {
  validatePosteriorRepresentation(representation, engine)
  stop(
    "posteriorDiagnostics is not yet implemented for engine '",
    posteriorEngineName(engine),
    "'",
    call. = FALSE
  )
}

#' Extract the representative posterior parameter value used by a fit.
#'
#' @param engine An internal `psPosteriorEngine` object.
#' @param model An internal `psModel` object.
#' @param representation A `psPosteriorRepresentation` object.
#' @param ... Engine- and model-specific controls.
#' @return Named posterior point estimate when implemented.
#' @keywords internal
#' @noRd
posteriorPointEstimate = function(engine, model, representation, ...) {
  UseMethod("posteriorPointEstimate")
}

#' @rdname posteriorPointEstimate
#' @keywords internal
#' @noRd
posteriorPointEstimate.psPosteriorEngine = function(engine,
                                                     model,
                                                     representation,
                                                     ...) {
  validateEngineModelPair(engine, model)
  validatePosteriorRepresentation(representation, engine)
  stop(
    "posteriorPointEstimate is not yet implemented for engine '",
    posteriorEngineName(engine),
    "'",
    call. = FALSE
  )
}

#' Validate that a model/engine combination is declared as supported.
#'
#' @param engine An internal `psPosteriorEngine` object.
#' @param model An internal `psModel` object.
#' @return The inputs, invisibly, when the pairing is supported.
#' @keywords internal
#' @noRd
validateEngineModelPair = function(engine, model) {
  if (!inherits(engine, "psPosteriorEngine")) {
    stop("engine must inherit from psPosteriorEngine")
  }
  if (!inherits(model, "psModel")) {
    stop("model must inherit from psModel")
  }

  if (!supportsPosteriorEngine(model, engine)) {
    stop(
      "Posterior engine '", posteriorEngineName(engine),
      "' is not supported by model '", model$model, "'",
      call. = FALSE
    )
  }

  invisible(list(engine = engine, model = model))
}
