#' S3 print method for an object of class \code{psFit}
#'
#' @param x an object of class \code{psFit}, usually from \code{\link{fitDist}}
#'   or \code{\link{fitZIDist}}.
#' @param nterms Number of probability terms to print. If \code{NULL}, print
#'   the terms already stored in the fitted object, capped at 10 for posterior
#'   and bootstrap summaries.
#' @param ... other arguments passed to delegated print methods.
#'
#' @return No return value, called for side effects.
#' @export
print.psFit = function(x, nterms = NULL, ...) {
  isBayes = identical(x$method, "bayes")

  if (isBayes && inherits(x$posterior, "psPosterior")) {
    print(x$posterior, nterms = nterms, ...)
    return(invisible(x))
  }

  if (x$model == "ziz") {
    cat("The estimated mixing parameter, pi, is", signif(x$pi, 4), "\n")
  }

  if (x$model %in% c("zeta", "ziz")) {
    cat("The estimated shape parameter is", round(x$shape, 4), "\n")
  }

  if (x$model %in% c("log", "logarithmic")) {
    cat("The estimated model parameter pi is", signif(x$pi, 4), "\n")
  }

  if (x$model == "zeta") {
    cat(
      "The standard error of shape parameter is",
      round(sqrt(x$var.shape), 4),
      "\n"
    )
  }

  if (x$model %in% c("zeta", "ziz", "log", "logarithmic")) {
    if (is.null(nterms)) {
      nterms = length(x$fitted)
    }
    if (!is.numeric(nterms) || length(nterms) != 1L || !is.finite(nterms) ||
        nterms < 1 || nterms != floor(nterms)) {
      stop("nterms must be one positive integer")
    }

    fittedValues = fitted(x, n = as.integer(nterms), type = "plugIn")
    cat("The first", nterms, "fitted values are:\n")
    print(fittedValues)
  }

  if (inherits(x$bootstrap, "psBootstrap")) {
    diagnostics = x$bootstrap$diagnostics
    cat("\nBootstrap uncertainty is attached")
    if (!is.null(diagnostics$B)) {
      cat(" (", diagnostics$B, " replicates)", sep = "")
    }
    cat(".\n")
    cat("Use bootstrapProbs(), summary(), or plot(x$bootstrap) for details.\n")
  }

  invisible(x)
}
