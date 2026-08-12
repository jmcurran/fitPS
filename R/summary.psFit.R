#' S3 summary method for an object of class \code{psFit}
#'
#' @param object an object of class \code{psFit}, usually from
#'   \code{\link{fitDist}} or \code{\link{fitZIDist}}.
#' @param nterms Optional number of leading posterior or bootstrap probability
#'   summaries to include when an uncertainty object is attached.
#' @param ... other arguments passed to delegated summary methods.
#'
#' @details
#' Bayesian fits delegate to the attached \code{psPosterior} object. For
#' maximum-likelihood fits, the existing parameter summary and zero-inflation
#' likelihood-ratio test are preserved. If a \code{psBootstrap} object is
#' attached, its parameter and probability summaries are printed afterwards.
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
    llH0 = -fitDist(object$psData)$fit$value
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
  } else {
    stop("summary is not currently implemented for this fitted model")
  }

  printCoefmat(cmat)

  if (inherits(object$bootstrap, "psBootstrap")) {
    cat("\n")
    print(summary(object$bootstrap, nterms = nterms, ...))
  }

  invisible(cmat)
}
