#' Signal use of a deprecated fitting interface
#'
#' Legacy fitting helpers remain callable for compatibility, but new code should
#' use the model-oriented `fit()` interface. The package option
#' `fitPS.deprecationWarnings` may be set to `FALSE` by automated test suites or
#' other controlled callers that deliberately exercise compatibility paths.
#'
#' @param old Name of the deprecated fitting helper.
#' @param replacement Replacement call to recommend.
#' @return Invisibly returns `NULL`.
#' @keywords internal
#' @noRd
signalDeprecatedInterface = function(old, replacement) {
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

#' Signal use of a deprecated fitting interface
#'
#' @param old Name of the deprecated fitting helper.
#' @param replacement Replacement call to recommend.
#' @return Invisibly returns `NULL`.
#' @keywords internal
#' @noRd
signalDeprecatedFitter = function(old, replacement) {
  signalDeprecatedInterface(old, replacement)
}
