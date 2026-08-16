#' Fit a fitPS model
#'
#' Fits a `psModel` object to forensic P- or S-survey data using maximum
#' likelihood, parametric Bayesian inference, the ordinary nonparametric
#' bootstrap, or Rubin's Bayesian Bootstrap. `fit()` is the canonical public
#' entry point for model fitting and uncertainty procedures.
#'
#' @param x An object of class `psData`.
#' @param model A model descriptor inheriting from `psModel`, such as an object
#'   returned by `zetaModel()`, `zizModel()`, or `logarithmicModel()`.
#' @param nterms Number of fitted P/S probability terms to retain.
#' @param method Inferential method. Use `"mle"` for maximum likelihood,
#'   `"bayes"` for parametric Bayesian inference, `"bootstrap"` for the
#'   observation-level nonparametric bootstrap, or `"bayesianBootstrap"` for
#'   Rubin's Bayesian Bootstrap. Legacy Bayesian engine names remain accepted
#'   and are normalised to `"bayes"`.
#' @param prior Optional prior used for parametric Bayesian fitting. For
#'   external Bayesian models this may be any model-specific object understood
#'   by `modelLogPrior()`.
#' @param bayesOptions Optional parametric Bayesian fitting controls.
#' @param B Number of bootstrap or Bayesian Bootstrap replicates.
#' @param level Confidence or equal-tail interval level used by bootstrap
#'   methods.
#' @param seed Optional random-number seed used by bootstrap methods.
#' @param silent Logical; suppress ordinary bootstrap progress messages when
#'   `TRUE`.
#' @param parallel Logical; use parallel ordinary-bootstrap fitting when
#'   `TRUE`.
#' @param progressBar Logical; display an ordinary-bootstrap progress bar when
#'   `TRUE`.
#' @param pbopts Options passed to `pbapply::pboptions()` when an ordinary
#'   bootstrap progress bar is requested.
#' @param ... Additional controls passed to the model-specific fitting path.
#' @return For `method = "mle"` or `"bayes"`, an object of class `psFit`.
#'   For `method = "bootstrap"`, a maximum-likelihood `psFit` with a
#'   `psBootstrap` object attached as `bootstrap`. For
#'   `method = "bayesianBootstrap"`, a `psBayesianBootstrap` object.
#' @references
#' \insertRef{curran2024}{fitPS}
#'
#' \insertRef{efron1979}{fitPS}
#'
#' \insertRef{rubin1981}{fitPS}
#'
#' @export
#'
#' @examples
#' data(Psurveys)
#' zetaFit = fit(Psurveys$roux, model = zetaModel())
#' logFit = fit(Psurveys$roux, model = logarithmicModel())
#' if (interactive()) {
#'   bootFit = fit(
#'     Psurveys$roux,
#'     model = zizModel(),
#'     method = "bootstrap",
#'     B = 20,
#'     seed = 123,
#'     silent = TRUE,
#'     parallel = FALSE
#'   )
#'   bayesBoot = fit(
#'     Psurveys$roux,
#'     model = zetaModel(),
#'     method = "bayesianBootstrap",
#'     B = 20,
#'     seed = 123
#'   )
#' }
fit = function(x,
               model,
               nterms = 10,
               method = c(
                 "mle", "bayes", "bootstrap", "bayesianBootstrap",
                 "integrate", "numerical", "mcmc", "laplace", "importance"
               ),
               prior,
               bayesOptions = NULL,
               B = 2000,
               level = 0.95,
               seed = NULL,
               silent = FALSE,
               parallel = TRUE,
               progressBar = FALSE,
               pbopts = list(type = "txt"),
               ...) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData")
  }
  if (!inherits(model, "psModel")) {
    stop("model must inherit from psModel")
  }

  method = match.arg(method)

  if (identical(method, "bootstrap")) {
    if (!model$model %in% c("zeta", "ziz")) {
      stop(
        "method = 'bootstrap' currently supports zeta and ziz models",
        call. = FALSE
      )
    }

    mleFit = do.call(
      fitModel,
      c(
        list(
          model = model,
          x = x,
          nterms = nterms,
          method = "mle"
        ),
        list(...)
      )
    )

    result = bootstrapPsFit(
      object = mleFit,
      B = B,
      level = level,
      seed = seed,
      silent = silent,
      parallel = parallel,
      progressBar = progressBar,
      pbopts = pbopts
    )
    result$method = "bootstrap"
    return(result)
  }

  if (identical(method, "bayesianBootstrap")) {
    return(bayesianBootstrapModel(
      x = x,
      model = model,
      B = B,
      seed = seed,
      nterms = nterms,
      level = level
    ))
  }

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
  if (!is.null(seed)) {
    args$seed = seed
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
