#' Compute zero-inflated zeta P or S probabilities
#'
#' Maps one or more pairs of zero-inflated zeta parameters to the requested
#' P- or S-survey probabilities. This is the shared internal transformation
#' used by posterior probability summaries.
#'
#' @param pi Numeric vector of mixing probabilities in the interval `[0, 1]`.
#' @param shape Numeric vector of standard zeta shape parameters greater than
#'   one.
#' @param n Non-negative integer P-term indices or positive integer S-term
#'   indices.
#' @param type Survey type, either `"P"` or `"S"`.
#'
#' @return A numeric matrix with one row per parameter pair and one column per
#'   requested probability. Columns are named with the corresponding P or S
#'   term.
#'
#' @keywords internal
zizProbabilities = function(pi, shape, n, type) {
  if (!is.numeric(pi) || length(pi) == 0L || any(!is.finite(pi)) ||
      any(pi < 0 | pi > 1)) {
    stop("pi must contain finite values in [0, 1]")
  }

  if (!is.numeric(shape) || length(shape) == 0L || any(!is.finite(shape))) {
    stop("shape must contain finite numeric values")
  }
  validateZetaShape(shape)

  parameterCount = max(length(pi), length(shape))
  if (!length(pi) %in% c(1L, parameterCount) ||
      !length(shape) %in% c(1L, parameterCount)) {
    stop("pi and shape must have equal lengths or one must have length one")
  }

  pi = rep(pi, length.out = parameterCount)
  shape = rep(shape, length.out = parameterCount)

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
  inflatedIndex = if (type == "P") 0L else 1L

  probabilities = vapply(
    seq_along(latentValues),
    function(termIndex) {
      values = (1 - pi) * dzetaStandard(
        latentValues[termIndex],
        shape = shape
      )
      if (n[termIndex] == inflatedIndex) {
        values = values + pi
      }
      values
    },
    numeric(parameterCount)
  )

  if (parameterCount == 1L) {
    probabilities = matrix(probabilities, nrow = 1L)
  }

  colnames(probabilities) = paste0(type, n)
  probabilities
}
