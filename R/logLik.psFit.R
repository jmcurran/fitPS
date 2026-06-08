#' Extract the log-likelihood from a zeta model fit
#'
#' Returns the maximised log-likelihood for a fitted zeta model.
#'
#' This method allows generic functions such as [stats::AIC()] and
#' [stats::BIC()] to work with objects of class `"psFit"`.
#'
#' @param object An object of class `"psFit"`.
#' @param ... Additional arguments passed to methods. Currently ignored.
#'
#' @return An object of class `"logLik"` with attributes `"df"` and `"nobs"`.
#'
#' @export
logLik.psFit = function(object, ...) {
  if (!is(object, "psFit")) {
    stop("object must be an object of class psFit")
  }

  if (is.null(object$fit$value)) {
    stop("object does not contain an optimised log-likelihood value")
  }

  value = -object$fit$value
  attr(value, "df") = length(object$fit$par)
  attr(value, "nobs") = sum(object$psData$data$rn)
  class(value) = "logLik"

  value
}
