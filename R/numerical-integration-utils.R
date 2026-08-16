#' Construct one-dimensional quadrature weights for an evaluated grid
#'
#' Prefer composite Simpson quadrature on an equally spaced grid. If an odd
#' final interval remains, use a trapezoidal contribution for that interval.
#' Non-equally-spaced grids fall back to composite trapezoidal weights.
#'
#' @param x Strictly increasing grid locations.
#' @return A list containing numeric quadrature weights and the rule used.
#' @keywords internal
#' @noRd
posteriorQuadratureWeights = function(x) {
  x = as.numeric(x)
  if (length(x) < 2L || any(!is.finite(x)) || any(diff(x) <= 0)) {
    stop("x must contain at least two finite, strictly increasing grid values", call. = FALSE)
  }

  differences = diff(x)
  spacing = differences[[1L]]
  tolerance = 100 * .Machine$double.eps * max(1, abs(spacing), max(abs(x)))
  equallySpaced = max(abs(differences - spacing)) <= tolerance

  if (!equallySpaced || length(x) < 3L) {
    weights = numeric(length(x))
    weights[[1L]] = differences[[1L]] / 2
    weights[[length(x)]] = differences[[length(differences)]] / 2
    if (length(x) > 2L) {
      weights[2L:(length(x) - 1L)] =
        (differences[-length(differences)] + differences[-1L]) / 2
    }
    return(list(weights = weights, method = "trapezoid"))
  }

  n = length(x)
  nIntervals = n - 1L
  nSimpsonIntervals = nIntervals - (nIntervals %% 2L)
  weights = numeric(n)

  if (nSimpsonIntervals >= 2L) {
    lastSimpsonPoint = nSimpsonIntervals + 1L
    coefficients = rep(2, lastSimpsonPoint)
    coefficients[c(1L, lastSimpsonPoint)] = 1
    coefficients[seq.int(2L, lastSimpsonPoint - 1L, by = 2L)] = 4
    weights[seq_len(lastSimpsonPoint)] = coefficients * spacing / 3
  }

  method = "simpson"
  if (nSimpsonIntervals < nIntervals) {
    lastInterval = x[[n]] - x[[n - 1L]]
    weights[[n - 1L]] = weights[[n - 1L]] + lastInterval / 2
    weights[[n]] = lastInterval / 2
    method = "simpson+trapezoid"
  }

  list(weights = weights, method = method)
}

#' Construct a cumulative integral over an evaluated one-dimensional grid
#'
#' For an equally spaced grid, cumulative areas are derived from the same local
#' quadratic interpolants used by Simpson's rule. If the grid has an unmatched
#' final interval, that final contribution is trapezoidal. Non-equally-spaced
#' grids use the composite trapezoidal rule.
#'
#' @param x Strictly increasing grid locations.
#' @param y Function values at `x`.
#' @return A list containing cumulative areas, total area, and integration rule.
#' @keywords internal
#' @noRd
cumulativePosteriorIntegral = function(x, y) {
  x = as.numeric(x)
  y = as.numeric(y)
  if (length(x) != length(y) || length(x) < 2L || any(!is.finite(x)) ||
      any(!is.finite(y)) || any(diff(x) <= 0)) {
    stop("x and y must be finite vectors of equal length on a strictly increasing grid", call. = FALSE)
  }

  quadrature = posteriorQuadratureWeights(x)
  cumulative = numeric(length(x))

  if (quadrature$method == "trapezoid") {
    increments = diff(x) * (y[-length(y)] + y[-1L]) / 2
    cumulative[-1L] = cumsum(increments)
  } else {
    n = length(x)
    nIntervals = n - 1L
    nSimpsonIntervals = nIntervals - (nIntervals %% 2L)

    for (startIndex in seq.int(1L, nSimpsonIntervals - 1L, by = 2L)) {
      h = x[[startIndex + 1L]] - x[[startIndex]]
      y0 = y[[startIndex]]
      y1 = y[[startIndex + 1L]]
      y2 = y[[startIndex + 2L]]

      firstHalf = h * (5 * y0 + 8 * y1 - y2) / 12
      pairArea = h * (y0 + 4 * y1 + y2) / 3
      cumulative[[startIndex + 1L]] = cumulative[[startIndex]] + firstHalf
      cumulative[[startIndex + 2L]] = cumulative[[startIndex]] + pairArea
    }

    if (nSimpsonIntervals < nIntervals) {
      lastIndex = n
      lastArea = (x[[lastIndex]] - x[[lastIndex - 1L]]) *
        (y[[lastIndex - 1L]] + y[[lastIndex]]) / 2
      cumulative[[lastIndex]] = cumulative[[lastIndex - 1L]] + lastArea
    }
  }

  area = sum(quadrature$weights * y)
  if (!is.finite(area) || area <= 0) {
    stop("The evaluated posterior density could not be normalized.", call. = FALSE)
  }

  # Simpson interpolation can show tiny floating-point non-monotonicity in a
  # far tail. Preserve a valid numerical CDF without changing material mass.
  cumulative = pmin(cummax(cumulative), area)
  cumulative[[length(cumulative)]] = area

  list(
    cumulative = cumulative,
    area = area,
    method = quadrature$method,
    weights = quadrature$weights
  )
}
