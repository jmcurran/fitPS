#' Compute a bootstrap distribution for a fitted fitPS model
#'
#' Compute and attach a nonparametric bootstrap distribution to a maximum
#' likelihood \code{psFit} object. The returned fitted object contains a
#' \code{psBootstrap} object in \code{object$bootstrap}.
#'
#' @param object A maximum likelihood \code{psFit} object.
#' @param B Number of bootstrap replicates.
#' @param level Confidence level used for percentile intervals.
#' @param seed Optional random seed for reproducible resampling.
#' @param silent Logical; suppress progress messages when \code{TRUE}.
#' @param parallel Logical; use parallel fitting when \code{TRUE}.
#' @param progressBar Logical; display a progress bar when \code{TRUE}.
#' @param pbopts Options passed to \code{pbapply::pboptions()} when progress
#'   bars are requested.
#' @param ... Additional arguments passed to methods.
#'
#' @return The input \code{psFit} object with a \code{psBootstrap} object
#'   attached as \code{bootstrap}.
#'
#' @details
#' The bootstrap is an observation-level nonparametric bootstrap. Each
#' successful parameter replicate is transformed into the corresponding P or S
#' probabilities before probability summaries are calculated. Consequently,
#' the reported bootstrap mean probabilities are averages of transformed
#' bootstrap replicates, not probabilities evaluated at average bootstrap
#' parameter values.
#'
#' @examples
#' if (interactive()) {
#' data(Psurveys)
#' fit = fitZIDist(Psurveys$roux, nterms = 4)
#' fit = bootstrapFit(
#'   fit,
#'   B = 20,
#'   seed = 123,
#'   silent = TRUE,
#'   parallel = FALSE
#' )
#' bootstrapProbs(fit)
#'
#' }
#' @export
bootstrapFit = function(object, ...) {
  UseMethod("bootstrapFit")
}

#' @rdname bootstrapFit
#' @export
bootstrapFit.psFit = function(object,
                               B = 2000,
                               level = 0.95,
                               seed = NULL,
                               silent = FALSE,
                               parallel = TRUE,
                               progressBar = FALSE,
                               pbopts = list(type = "txt"),
                               ...) {
  bootstrapPsFit(
    object = object,
    B = B,
    level = level,
    seed = seed,
    silent = silent,
    parallel = parallel,
    progressBar = progressBar,
    pbopts = pbopts
  )
}
