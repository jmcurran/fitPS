#' Convert P/S data to the positive-integer support used by the ZIZ likelihood.
#'
#' @param x An input object or numeric vector required by the helper.
#' @return A numeric vector on the positive-integer ZIZ support.
#' @keywords internal
#' @noRd
zizObservationData = function(x) {
  psObservationData(x)
}

#' Evaluate the zero-inflated zeta log likelihood for aggregated observations.
#'
#' @param obsData Positive-integer observation support used by the ZIZ likelihood.
#' @param counts Observation frequencies corresponding to `obsData`.
#' @param pi Zero-inflation probability on the natural scale.
#' @param shape Zeta shape parameter on the fitPS scale.
#' @return A numeric log-likelihood value, or `-Inf` outside the parameter space.
#' @keywords internal
#' @noRd
zizLogLikelihood = function(obsData, counts, pi, shape) {
  pi = unname(pi)
  shape = unname(shape)

  if (!is.numeric(pi) || length(pi) != 1L || !is.finite(pi) || pi <= 0 || pi >= 1) {
    return(-Inf)
  }

  if (!is.numeric(shape) || length(shape) != 1L || !is.finite(shape) || shape <= 1) {
    return(-Inf)
  }

  probabilities = (1 - pi) * dzetaStandard(obsData, shape = shape)
  probabilities[obsData == 1] = probabilities[obsData == 1] + pi

  if (any(!is.finite(probabilities)) || any(probabilities <= 0)) {
    return(-Inf)
  }

  sum(counts * log(probabilities))
}

#' Construct and normalise a two-dimensional numerical posterior grid for ZIZ parameters.
#'
#' @param x An input object or numeric vector required by the helper.
#' @param prior A prior object created by `makePrior()`.
#' @param shape1 First beta-prior shape parameter for the zero-inflation probability.
#' @param shape2 Second beta-prior shape parameter for the zero-inflation probability.
#' @param nPiGrid Number of grid points for the zero-inflation probability.
#' @param nShapeGrid Number of grid points for the zeta shape parameter.
#' @return A list describing the normalised joint grid, marginals, moments, and grid spacings.
#' @keywords internal
#' @noRd
makeZizPosteriorGrid = function(x,
                                prior,
                                shape1 = 1,
                                shape2 = 1,
                                nPiGrid = 101,
                                nShapeGrid = 101) {
  validateBayesPrior(prior)

  if (!is.numeric(shape1) || length(shape1) != 1L || !is.finite(shape1) || shape1 <= 0) {
    stop("shape1 must be a positive finite number")
  }

  if (!is.numeric(shape2) || length(shape2) != 1L || !is.finite(shape2) || shape2 <= 0) {
    stop("shape2 must be a positive finite number")
  }

  nPiGrid = as.integer(nPiGrid)
  nShapeGrid = as.integer(nShapeGrid)

  if (!is.finite(nPiGrid) || nPiGrid < 9L) {
    stop("nPiGrid must be at least 9")
  }

  if (!is.finite(nShapeGrid) || nShapeGrid < 9L) {
    stop("nShapeGrid must be at least 9")
  }

  obsData = zizObservationData(x)
  counts = x$data$rn
  # Keep the numerical grid just inside the open support of pi. The clipping
  # avoids log-density failures at exactly 0 or 1 without changing the model.
  piEps = sqrt(.Machine$double.eps)
  piGrid = seq(piEps, 1 - piEps, length.out = nPiGrid)
  shapeGrid = seq(prior$range[1], prior$range[2], length.out = nShapeGrid)

  logPosterior = outer(
    piGrid,
    shapeGrid,
    Vectorize(function(pi, shape) {
      zizLogLikelihood(obsData, counts, pi, shape) +
        dbeta(pi, shape1, shape2, log = TRUE) +
        prior$logd(shape)
    })
  )

  finiteValues = logPosterior[is.finite(logPosterior)]
  if (length(finiteValues) == 0L) {
    stop("The zero-inflated zeta posterior grid has no finite posterior values")
  }

  # Subtract the maximum log posterior before exponentiating. This standard
  # log-sum-exp rescaling preserves relative weights while avoiding overflow
  # and severe underflow during numerical integration.
  logScale = max(finiteValues)
  scaledWeights = exp(logPosterior - logScale)
  scaledWeights[!is.finite(scaledWeights)] = 0

  # The grid is rectangular and equally spaced, so each cell has integration
  # weight dPi * dShape. Marginal densities below retain the complementary
  # spacing factor so subsequent one-dimensional sums approximate integrals.
  dPi = diff(range(piGrid)) / (length(piGrid) - 1L)
  dShape = diff(range(shapeGrid)) / (length(shapeGrid) - 1L)
  scaledIntegral = sum(scaledWeights) * dPi * dShape

  if (!is.finite(scaledIntegral) || scaledIntegral <= 0) {
    stop("The zero-inflated zeta posterior grid could not be normalized")
  }

  jointDensity = scaledWeights / scaledIntegral
  marginalPiDensity = rowSums(jointDensity) * dShape
  marginalShapeDensity = colSums(jointDensity) * dPi

  piMean = sum(piGrid * marginalPiDensity) * dPi
  shapeMean = sum(shapeGrid * marginalShapeDensity) * dShape
  piVariance = sum((piGrid - piMean)^2 * marginalPiDensity) * dPi
  shapeVariance = sum((shapeGrid - shapeMean)^2 * marginalShapeDensity) * dShape
  covariance = sum(
    outer(piGrid - piMean, shapeGrid - shapeMean) * jointDensity
  ) * dPi * dShape

  list(
    pi = piGrid,
    shape = shapeGrid,
    logPosterior = logPosterior,
    density = jointDensity,
    marginalDensity = list(
      pi = marginalPiDensity,
      shape = marginalShapeDensity
    ),
    normalizingConstant = exp(logScale) * scaledIntegral,
    dPi = dPi,
    dShape = dShape,
    mean = c(pi = piMean, shape = shapeMean),
    varCov = matrix(
      c(piVariance, covariance, covariance, shapeVariance),
      nrow = 2L,
      dimnames = list(c("pi", "shape"), c("pi", "shape"))
    )
  )
}
