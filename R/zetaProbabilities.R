#' Compute zeta P or S probabilities
#'
#' Maps one or more standard zeta shape parameters to the requested P- or
#' S-survey probabilities. This is the frequentist bootstrap counterpart to
#' the zero-inflated transformation in `zizProbabilities()`.
#'
#' @param shape Numeric vector of standard zeta shape parameters greater than
#'   one.
#' @param n Non-negative integer P-term indices or positive integer S-term
#'   indices.
#' @param type Survey type, either `"P"` or `"S"`.
#'
#' @return A numeric matrix with one row per shape value and one column per
#'   requested probability. Columns are named with the corresponding P or S
#'   term.
#'
#' @keywords internal
zetaProbabilities = function(shape, n, type) {
  if (!is.numeric(shape) || length(shape) == 0L || any(!is.finite(shape))) {
    stop("shape must contain finite numeric values")
  }
  validateZetaShape(shape)

  if (!is.numeric(n) || length(n) == 0L || any(!is.finite(n)) ||
      any(n != floor(n))) {
    stop("n must contain finite integer values")
  }
  n = as.integer(n)

  type = match.arg(type, c("P", "S"))
  if (type == "P" && any(n < 0L)) {
    stop("P probability indices must be non-negative")
  }
  if (type == "S" && any(n < 1L)) {
    stop("S probability indices must be positive")
  }

  latentValues = if (type == "P") {
    n + 1L
  } else {
    n
  }

  probabilities = vapply(
    latentValues,
    function(value) {
      dzetaStandard(value, shape = shape)
    },
    numeric(length(shape))
  )

  if (length(shape) == 1L) {
    probabilities = matrix(probabilities, nrow = 1L)
  }

  colnames(probabilities) = paste0(type, n)
  probabilities
}
