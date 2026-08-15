#' Fit a fitPS model
#'
#' Fits a `psModel` object to forensic P- or S-survey data. This is the common
#' model-oriented fitting entry point introduced for the fitPS public model
#' API. Stage 8.2 supports the built-in zeta, zero-inflated zeta, and
#' logarithmic model descriptors. The public third-party model contract is
#' introduced separately so it can be validated against an external model.
#'
#' @param x An object of class `psData`.
#' @param model A model descriptor inheriting from `psModel`, such as an object
#'   returned by `zetaModel()`, `zizModel()`, or `logarithmicModel()`.
#' @param nterms Number of fitted P/S probability terms to retain.
#' @param method Fitting method. Built-in models currently preserve the methods
#'   accepted by their established fitting functions.
#' @param prior Optional prior used for Bayesian fitting.
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

#' Fit one model descriptor through its established fitting implementation.
#'
#' Stage 8.2 keeps this dispatcher internal while the third-party MLE contract
#' is validated. Built-in methods deliberately delegate to the established
#' fitters so `fit()` does not duplicate their fitting mathematics.
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
#' @exportS3Method fitModel zetaModel
#' @noRd
fitModel.zetaModel = function(model, x, ...) {
  result = fitDist(x = x, ...)
  result$modelObject = model
  result
}

#' @rdname fitModel
#' @keywords internal
#' @exportS3Method fitModel zizModel
#' @noRd
fitModel.zizModel = function(model, x, ...) {
  result = fitZIDist(x = x, ...)
  result$modelObject = model
  result
}

#' @rdname fitModel
#' @keywords internal
#' @exportS3Method fitModel logarithmicModel
#' @noRd
fitModel.logarithmicModel = function(model, x, ...) {
  result = fitlogDist(x = x, ...)
  result$modelObject = model
  result
}

#' @rdname fitModel
#' @keywords internal
#' @exportS3Method fitModel psModel
#' @noRd
fitModel.psModel = function(model, x, ...) {
  stop(
    "fitting is not yet implemented for model class '",
    class(model)[1L],
    "'; the third-party model fitting contract is introduced in Stage 8.3",
    call. = FALSE
  )
}
