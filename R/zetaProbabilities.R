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

  type = normaliseSurveyType(type)
  n = normaliseProbabilityIndices(n, type)
  latentValues = latentPsValues(n, type)

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

  colnames(probabilities) = psProbabilityTermNames(n, type)
  probabilities
}
