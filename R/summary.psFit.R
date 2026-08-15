#' S3 summary method for an object of class \code{psFit}
#'
#' @param object An object of class \code{psFit}, usually returned by [fit()] or
#'   one of the deprecated compatibility fitters.
#' @param nterms Optional number of leading posterior or bootstrap probability
#'   summaries to include when an uncertainty object is attached.
#' @param ... Other arguments passed to delegated summary methods.
#'
#' @details
#' Bayesian fits delegate to the attached \code{psPosterior} object. For
#' maximum-likelihood fits, built-in zeta, ZIZ, and logarithmic presentation is
#' preserved. Other models use the public model contract to report their fitted
#' natural-scale parameters and, when available, standard errors obtained from
#' the fitted covariance matrix. If a \code{psBootstrap} object is attached,
#' its parameter and probability summaries are printed afterwards.
#'
#' @importFrom stats pchisq printCoefmat
#' @return For Bayesian fits, a \code{summary.psPosterior} object. For
#'   maximum-likelihood fits, invisibly returns the coefficient matrix.
#' @export
summary.psFit = function(object, nterms = NULL, ...) {
  if (identical(object$method, "bayes") &&
      inherits(object$posterior, "psPosterior")) {
    return(summary(object$posterior, nterms = nterms, ...))
  }

  if (object$model == "zeta") {
    cmat = matrix(c(object$shape, sqrt(object$var.shape)), nrow = 1)
    colnames(cmat) = c("Estimate", "Std.Err")
    rownames(cmat) = "shape"
  } else if (object$model == "ziz") {
    llH0 = -fit(object$psData, model = zetaModel())$fit$value
    llMle = -object$fit$value
    lrtStat = 2 * (llMle - llH0)

    cmat = rbind(
      c(object$shape, sqrt(object$var.cov[2, 2]), NA, NA),
      c(
        object$pi,
        sqrt(object$var.cov[1, 1]),
        lrtStat,
        pchisq(lrtStat, 1, lower.tail = FALSE)
      )
    )
    colnames(cmat) = c("Estimate", "Std.Err", "X^2 value", "Pr(>X^2)")
    rownames(cmat) = c("shape", "pi")
  } else if (object$model %in% c("log", "logarithmic")) {
    cmat = matrix(c(object$pi, sqrt(object$var.pi)), nrow = 1L)
    colnames(cmat) = c("Estimate", "Std.Err")
    rownames(cmat) = "pi"
  } else {
    model = modelFromFit(object)
    estimates = fitModelParameters(object, model)
    standardErrors = rep(NA_real_, length(estimates))
    names(standardErrors) = names(estimates)

    covariance = object$var.cov
    if (is.matrix(covariance) &&
        all(names(estimates) %in% rownames(covariance)) &&
        all(names(estimates) %in% colnames(covariance))) {
      variances = diag(covariance[names(estimates), names(estimates), drop = FALSE])
      valid = is.finite(variances) & variances >= 0
      standardErrors[valid] = sqrt(variances[valid])
    }

    cmat = cbind(
      Estimate = as.numeric(estimates),
      Std.Err = as.numeric(standardErrors)
    )
    rownames(cmat) = names(estimates)
  }

  printCoefmat(cmat)

  if (inherits(object$bootstrap, "psBootstrap")) {
    cat("\n")
    print(summary(object$bootstrap, nterms = nterms, ...))
  }

  invisible(cmat)
}
