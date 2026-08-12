#' Summarise posterior zero-inflated zeta probabilities
#'
#' Converts posterior parameter representations into a common set of posterior
#' summaries for P or S probabilities.
#'
#' @param probabilities Numeric matrix with posterior realisations in rows and
#'   probability terms in columns.
#' @param weights Optional non-negative posterior weights.
#' @param level Equal-tailed credible interval level.
#' @param posteriorMethod Posterior engine used to produce the representation.
#'
#' @return A data frame containing posterior means, standard deviations, and
#'   equal-tailed credible intervals.
#'
#' @keywords internal
summariseZizProbabilities = function(probabilities,
                                     weights = NULL,
                                     level = 0.95,
                                     posteriorMethod = NA_character_) {
  probabilities = as.matrix(probabilities)

  if (!is.numeric(probabilities) || nrow(probabilities) == 0L ||
      ncol(probabilities) == 0L || any(!is.finite(probabilities))) {
    stop("probabilities must be a non-empty finite numeric matrix")
  }

  if (is.null(colnames(probabilities)) || any(colnames(probabilities) == "")) {
    stop("probabilities must have non-empty column names")
  }

  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("level must be one finite number strictly between zero and one")
  }

  if (is.null(weights)) {
    weights = rep(1 / nrow(probabilities), nrow(probabilities))
  } else {
    if (!is.numeric(weights) || length(weights) != nrow(probabilities) ||
        any(!is.finite(weights)) || any(weights < 0)) {
      stop("weights must be finite, non-negative, and match the number of rows")
    }

    weightSum = sum(weights)
    if (!is.finite(weightSum) || weightSum <= 0) {
      stop("weights must have a positive finite sum")
    }
    weights = weights / weightSum
  }

  lowerProbability = (1 - level) / 2
  upperProbability = 1 - lowerProbability
  estimates = colSums(probabilities * weights)
  secondMoments = colSums(probabilities^2 * weights)
  variances = pmax(0, secondMoments - estimates^2)

  intervals = vapply(
    seq_len(ncol(probabilities)),
    function(columnIndex) {
      weightedZizQuantile(
        values = probabilities[, columnIndex],
        weights = weights,
        probabilities = c(lowerProbability, upperProbability)
      )
    },
    numeric(2L)
  )

  data.frame(
    term = colnames(probabilities),
    estimate = unname(estimates),
    sd = sqrt(unname(variances)),
    lower = unname(intervals[1L, ]),
    upper = unname(intervals[2L, ]),
    level = rep(level, ncol(probabilities)),
    posteriorMethod = rep(posteriorMethod, ncol(probabilities)),
    stringsAsFactors = FALSE
  )
}

#' Weighted quantiles for posterior probability summaries
#'
#' @param values Numeric values.
#' @param weights Non-negative weights.
#' @param probabilities Quantile probabilities in `[0, 1]`.
#'
#' @return Numeric weighted quantiles.
#'
#' @keywords internal
weightedZizQuantile = function(values, weights, probabilities) {
  if (!is.numeric(values) || !is.numeric(weights) ||
      length(values) != length(weights) || length(values) == 0L ||
      any(!is.finite(values)) || any(!is.finite(weights)) || any(weights < 0)) {
    stop("values and weights must be finite numeric vectors of equal length")
  }

  if (!is.numeric(probabilities) || any(!is.finite(probabilities)) ||
      any(probabilities < 0 | probabilities > 1)) {
    stop("probabilities must lie in [0, 1]")
  }

  weightSum = sum(weights)
  if (!is.finite(weightSum) || weightSum <= 0) {
    stop("weights must have a positive finite sum")
  }

  orderIndex = order(values)
  sortedValues = values[orderIndex]
  sortedWeights = weights[orderIndex] / weightSum
  cumulativeWeights = cumsum(sortedWeights)

  vapply(
    probabilities,
    function(probability) {
      if (probability <= 0) {
        return(sortedValues[1L])
      }
      if (probability >= 1) {
        return(sortedValues[length(sortedValues)])
      }

      quantileIndex = which(cumulativeWeights >= probability)[1L]
      sortedValues[quantileIndex]
    },
    numeric(1L)
  )
}

#' Return the support indices used to label fitted P or S probabilities.
#'
#' @param type P- or S-survey type.
#' @param nterms Number of fitted P/S probability terms to retain.
#' @return Integer support indices for naming P/S probability terms.
#' @keywords internal
#' @noRd
zizProbabilityIndices = function(type, nterms) {
  type = match.arg(type, c("P", "S"))
  nterms = as.integer(nterms)

  if (!is.finite(nterms) || length(nterms) != 1L || nterms <= 0L) {
    stop("nterms must be a positive integer")
  }

  if (type == "P") {
    0:(nterms - 1L)
  } else {
    seq_len(nterms)
  }
}

#' Summarise fitted P/S probabilities over a numerical ZIZ posterior grid.
#'
#' @param posteriorGrid Numerical posterior-grid representation.
#' @param type P- or S-survey type.
#' @param nterms Number of fitted P/S probability terms to retain.
#' @param level Probability level for intervals or summaries.
#' @return A data frame of posterior fitted-probability summaries.
#' @keywords internal
#' @noRd
summariseZizGridProbabilities = function(posteriorGrid,
                                          type,
                                          nterms,
                                          level = 0.95) {
  piValues = rep(posteriorGrid$pi, times = length(posteriorGrid$shape))
  shapeValues = rep(posteriorGrid$shape, each = length(posteriorGrid$pi))
  weights = as.vector(posteriorGrid$density) *
    posteriorGrid$dPi * posteriorGrid$dShape

  probabilities = zizProbabilities(
    pi = piValues,
    shape = shapeValues,
    n = zizProbabilityIndices(type, nterms),
    type = type
  )

  summariseZizProbabilities(
    probabilities = probabilities,
    weights = weights,
    level = level,
    posteriorMethod = "numerical"
  )
}

#' Summarise fitted P/S probabilities over sampled ZIZ posterior representations.
#'
#' @param pi Zero-inflation probability on the natural scale.
#' @param shape Zeta shape parameter on the fitPS scale.
#' @param type P- or S-survey type.
#' @param nterms Number of fitted P/S probability terms to retain.
#' @param weights Non-negative normalised or normalisable weights.
#' @param level Probability level for intervals or summaries.
#' @param posteriorMethod Character label identifying the posterior approximation method.
#' @return A data frame of posterior fitted-probability summaries.
#' @keywords internal
#' @noRd
summariseZizSampleProbabilities = function(pi,
                                            shape,
                                            type,
                                            nterms,
                                            weights = NULL,
                                            level = 0.95,
                                            posteriorMethod) {
  probabilities = zizProbabilities(
    pi = pi,
    shape = shape,
    n = zizProbabilityIndices(type, nterms),
    type = type
  )

  summariseZizProbabilities(
    probabilities = probabilities,
    weights = weights,
    level = level,
    posteriorMethod = posteriorMethod
  )
}
