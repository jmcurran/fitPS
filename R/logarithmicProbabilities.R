#' Evaluate logarithmic P/S probabilities
#'
#' This helper evaluates the logarithmic-series distribution on the latent
#' positive-integer support used by fitPS. P probabilities use `n + 1`, while
#' S probabilities use `n` directly.
#'
#' @param pi Numeric logarithmic-distribution parameter values in `(0, 1)`.
#' @param n Requested P/S probability indices.
#' @param type Survey type, either `"P"` or `"S"`.
#' @return Numeric matrix with one row per parameter value and one column per
#'   requested probability term.
#' @importFrom VGAM dlog
#' @keywords internal
#' @noRd
logarithmicProbabilities = function(pi, n, type) {
  if (!is.numeric(pi) || length(pi) == 0L || any(!is.finite(pi))) {
    stop("pi must contain finite numeric values")
  }
  if (any(pi <= 0 | pi >= 1)) {
    stop("pi must lie strictly between 0 and 1")
  }

  type = normaliseSurveyType(type)
  n = normaliseProbabilityIndices(n, type)
  support = latentPsValues(n, type)

  values = outer(
    pi,
    support,
    Vectorize(function(parameter, value) {
      dlog(value, shape = parameter)
    })
  )
  colnames(values) = psProbabilityTermNames(n, type)
  values
}

#' Validate the logarithmic-distribution parameter
#'
#' @param pi Candidate scalar parameter.
#' @param name Name used in error messages.
#' @return `pi`, invisibly, when valid.
#' @keywords internal
#' @noRd
validateLogarithmicPi = function(pi, name = "pi") {
  if (!is.numeric(pi) || length(pi) != 1L || !is.finite(pi) ||
      pi <= 0 || pi >= 1) {
    stop(name, " must be one finite numeric value strictly between 0 and 1")
  }
  invisible(pi)
}
