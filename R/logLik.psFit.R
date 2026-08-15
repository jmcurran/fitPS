#' Extract the maximised log-likelihood from a fitPS model fit
#'
#' Returns the model log-likelihood evaluated at the maximum-likelihood estimate.
#' The returned `logLik` object includes the model parameter count and number of
#' observations, allowing [stats::AIC()] and [stats::BIC()] to use the shared
#' fitPS model-comparison contract.
#'
#' Bayesian fits are rejected because their stored representative parameter
#' values are posterior summaries rather than maximum-likelihood estimates.
#'
#' @param object An object of class `psFit`.
#' @param ... Additional arguments retained for S3 compatibility.
#' @return An object of class `logLik` with `df` and `nobs` attributes.
#' @export
logLik.psFit = function(object, ...) {
  if (!is(object, "psFit")) {
    stop("object must be an object of class psFit")
  }
  if (!identical(object$method, "mle")) {
    stop("AIC and BIC require a maximum-likelihood psFit object", call. = FALSE)
  }

  model = modelFromFit(object)
  parameters = fitModelParameters(object, model)
  value = modelLogLikelihood(model, parameters = parameters, data = object$psData)
  attr(value, "df") = modelParameterCount(model)
  attr(value, "nobs") = sum(object$psData$data$rn)
  class(value) = "logLik"
  value
}
