#' Signal use of a deprecated distribution-specific fitter
#'
#' The legacy fitters remain callable for compatibility, but new code should
#' use the model-oriented `fit()` interface. The package option
#' `fitPS.deprecationWarnings` may be set to `FALSE` by automated test suites or
#' other controlled callers that deliberately exercise compatibility paths.
#'
#' @param old Name of the deprecated fitter.
#' @param replacement Replacement call to recommend.
#' @return Invisibly returns `NULL`.
#' @keywords internal
#' @noRd
signalDeprecatedFitter = function(old, replacement) {
  if (!isTRUE(getOption("fitPS.deprecationWarnings", TRUE))) {
    return(invisible(NULL))
  }

  warning(
    old,
    "() is deprecated; use ",
    replacement,
    " instead.",
    call. = FALSE
  )
  invisible(NULL)
}
