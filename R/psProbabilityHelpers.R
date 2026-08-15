#' Validate and normalise P/S probability indices.
#'
#' @param n Candidate P- or S-term indices.
#' @param type Survey type, either `"P"` or `"S"`.
#' @return Integer probability indices.
#' @keywords internal
#' @noRd
normaliseProbabilityIndices = function(n, type) {
  if (!is.numeric(n) || length(n) == 0L || any(!is.finite(n)) ||
      any(n != floor(n))) {
    stop("n must contain finite integer values")
  }

  n = as.integer(n)
  type = normaliseSurveyType(type)

  if (type == "P" && any(n < 0L)) {
    stop("P probability indices must be non-negative")
  }
  if (type == "S" && any(n < 1L)) {
    stop("S probability indices must be positive")
  }

  n
}

#' Normalise a fitPS survey type.
#'
#' @param type Survey type, either `"P"` or `"S"`.
#' @return Normalised survey type.
#' @keywords internal
#' @noRd
normaliseSurveyType = function(type) {
  match.arg(type, c("P", "S"))
}

#' Convert P/S term indices to the positive-integer latent support.
#'
#' @param n Validated or candidate P/S probability indices.
#' @param type Survey type, either `"P"` or `"S"`.
#' @return Positive-integer support values used by the zeta-family likelihoods.
#' @keywords internal
#' @noRd
latentPsValues = function(n, type) {
  type = normaliseSurveyType(type)
  n = normaliseProbabilityIndices(n, type)

  if (type == "P") {
    n + 1L
  } else {
    n
  }
}

#' Construct standard fitPS probability term names.
#'
#' @param n Validated or candidate P/S probability indices.
#' @param type Survey type, either `"P"` or `"S"`.
#' @return Character vector such as `P0`, `P1` or `S1`, `S2`.
#' @keywords internal
#' @noRd
psProbabilityTermNames = function(n, type) {
  type = normaliseSurveyType(type)
  n = normaliseProbabilityIndices(n, type)
  paste0(type, n)
}

#' Convert aggregated fitPS observations to latent positive-integer support.
#'
#' @param x An object of class `psData`.
#' @return Integer observation support used by zeta-family likelihoods.
#' @keywords internal
#' @noRd
psObservationData = function(x) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData")
  }

  latentPsValues(x$data$n, x$type)
}

#' Validate the historical MLE support requirement.
#'
#' @param x An object of class `psData`.
#' @return `TRUE` invisibly when the MLE support requirement is satisfied.
#' @keywords internal
#' @noRd
validateMleObservationSupport = function(x) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData")
  }

  if (length(x$data$n) < 2L) {
    if (x$type == "S") {
      stop("There has to be at least one value higher than 1")
    } else {
      stop("There has to be at least one value higher than 0")
    }
  }

  invisible(TRUE)
}

#' Warn when Bayesian fitting has only one occupied support value.
#'
#' A proper prior may still produce a proper posterior when the likelihood has
#' no finite interior maximum. The warning makes the resulting prior sensitivity
#' explicit without treating sparse support as invalid data.
#'
#' @param x An object of class `psData`.
#' @return `NULL` invisibly.
#' @keywords internal
#' @noRd
warnSparseBayesianSupport = function(x) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData")
  }

  if (length(x$data$n) < 2L) {
    warning(
      "Survey data contain only one occupied support value. Bayesian fitting ",
      "is proceeding, but posterior inference may be strongly influenced by ",
      "the prior.",
      call. = FALSE
    )
  }

  invisible(NULL)
}
