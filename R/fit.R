#' Fit a fitPS model
#'
#' Fits a `psModel` object to forensic P- or S-survey data. Built-in models
#' retain their established fitting implementations. External models may use
#' the public model contract for fitPS-owned maximum-likelihood optimisation and
#' generic numerical or MCMC Bayesian fitting without modifying or rebuilding fitPS.
#'
#' @param x An object of class `psData`.
#' @param model A model descriptor inheriting from `psModel`, such as an object
#'   returned by `zetaModel()`, `zizModel()`, or `logarithmicModel()`.
#' @param nterms Number of fitted P/S probability terms to retain.
#' @param method Fitting method. External models support `"mle"` and generic
#'   `"bayes"` fitting when they advertise numerical or MCMC posterior engines.
#'   Legacy `"mcmc"` selection remains accepted and is normalised to Bayesian fitting.
#' @param prior Optional prior used for Bayesian fitting. For external Bayesian
#'   models this may be any model-specific object understood by `modelLogPrior()`.
#' @param bayesOptions Optional Bayesian fitting controls.
#' @param ... Additional controls passed to the model-specific fitting path.
#' @return An object of class `psFit` retaining the originating model descriptor
#'   in `modelObject` and the established character identifier in `model`.
#' @export
#'
#' @examples
#' data(Psurveys)
#' zetaFit = fit(Psurveys$roux, model = zetaModel())
#' logFit = fit(Psurveys$roux, model = logarithmicModel())
fit = function(x,
               model,
               nterms = 10,
               method = c(
                 "mle", "bayes", "integrate", "numerical", "mcmc",
                 "laplace", "importance"
               ),
               prior,
               bayesOptions = NULL,
               ...) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData")
  }
  if (!inherits(model, "psModel")) {
    stop("model must inherit from psModel")
  }

  method = match.arg(method)

  args = list(
    model = model,
    x = x,
    nterms = nterms,
    method = method,
    bayesOptions = bayesOptions,
    ...
  )
  if (!missing(prior)) {
    args$prior = prior
  }

  do.call(fitModel, args)
}

#' Select the default Bayesian engine for an external model.
#'
#' Deterministic numerical fitting is preferred for one- and two-parameter
#' models that advertise it. Models with three or more parameters use MCMC.
#'
#' @param model A `psModel` descriptor.
#' @return One posterior-engine name.
#' @keywords internal
#' @noRd
defaultExternalPosteriorMethod = function(model) {
  parameterCount = length(modelParameterNames(model))
  if (parameterCount <= 2L && supportsPosteriorEngine(model, "numerical")) {
    return("numerical")
  }
  if (supportsPosteriorEngine(model, "mcmc")) {
    return("mcmc")
  }

  supported = supportedPosteriorEngines(model)
  if (length(supported) == 0L) {
    stop("model does not declare a supported Bayesian posterior engine", call. = FALSE)
  }
  if (parameterCount > 2L && "numerical" %in% supported) {
    stop(
      "Models with three or more parameters require MCMC for Bayesian fitting",
      call. = FALSE
    )
  }
  supported[[1L]]
}

#' Fit one model descriptor through its established fitting implementation.
#'
#' Built-in methods delegate to the established internal fitting implementations
#' shared by `fit()` and the deprecated compatibility wrappers. The `psModel` fallback
#' provides the public generic MLE path for externally defined models.
#'
#' @param model A `psModel` descriptor.
#' @param x An object of class `psData`.
#' @param ... Arguments passed to the model-specific fitting implementation.
#' @return An object of class `psFit`.
#' @keywords internal
#' @noRd
fitModel = function(model, x, ...) {
  UseMethod("fitModel")
}
#' @rdname fitModel
#' @keywords internal
#' @exportS3Method fitModel psModel
#' @noRd
fitModel.psModel = function(model,
                             x,
                             nterms = 10,
                             method = "mle",
                             prior,
                             bayesOptions = NULL,
                             ...) {
  methodInfo = normaliseBayesMethod(method, bayesOptions = bayesOptions)
  method = methodInfo$method
  bayesOptions = methodInfo$bayesOptions

  if (identical(method, "bayes")) {
    automaticPosteriorMethod = defaultExternalPosteriorMethod(model)

    options = if (missing(prior)) {
      normaliseExternalBayesOptions(
        bayesOptions = bayesOptions,
        defaultPosteriorMethod = automaticPosteriorMethod
      )
    } else {
      normaliseExternalBayesOptions(
        bayesOptions = bayesOptions,
        prior = prior,
        defaultPosteriorMethod = automaticPosteriorMethod
      )
    }

    engine = posteriorEngine(options$posteriorMethod)
    validateEngineModelPair(engine, model)

    result = fitBayesianModel(
      model = model,
      posteriorMethod = options$posteriorMethod,
      x = x,
      prior = options$prior,
      nterms = nterms,
      ...
    )
    result$bayesOptions = options
    return(result)
  }

  validateMleObservationSupport(x)
  control = modelMleControl(model, x, ...)
  parameterNames = modelParameterNames(model)
  objective = function(parameterValues) {
    names(parameterValues) = parameterNames
    value = modelLogLikelihood(
      model = model,
      parameters = parameterValues,
      data = x
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
      "generic MLE optimisation failed for model '",
      model$model,
      "': ",
      fittedModel$message,
      call. = FALSE
    )
  }

  estimates = as.numeric(fittedModel$par)
  names(estimates) = parameterNames
  probabilityIndices = if (identical(x$type, "P")) {
    seq.int(0L, as.integer(nterms) - 1L)
  } else {
    seq_len(as.integer(nterms))
  }
  fittedMatrix = modelProbabilities(
    model = model,
    parameters = as.list(estimates),
    n = probabilityIndices,
    type = x$type
  )
  fittedValues = as.numeric(fittedMatrix[1L, ])
  names(fittedValues) = colnames(fittedMatrix)

  result = list(
    psData = x,
    fit = fittedModel,
    fitted = fittedValues,
    model = model$model,
    modelObject = model,
    method = "mle"
  )
  for (parameterName in parameterNames) {
    result[[parameterName]] = estimates[[parameterName]]
  }

  if (is.matrix(fittedModel$hessian) &&
      nrow(fittedModel$hessian) == length(parameterNames)) {
    covariance = tryCatch(
      solve(fittedModel$hessian),
      error = function(error) NULL
    )
    if (!is.null(covariance)) {
      rownames(covariance) = parameterNames
      colnames(covariance) = parameterNames
      result$var.cov = covariance
    }
  }

  class(result) = "psFit"
  result
}
