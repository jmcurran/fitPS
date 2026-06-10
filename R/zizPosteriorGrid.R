zizObservationData = function(x) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData")
  }

  if (length(x$data$n) < 2) {
    if (x$type == "S") {
      stop("There has to be at least one value higher than 1")
    } else {
      stop("There has to be at least one value higher than 0")
    }
  }

  if (x$type == "P") {
    x$data$n + 1
  } else {
    x$data$n
  }
}

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

  logScale = max(finiteValues)
  scaledWeights = exp(logPosterior - logScale)
  scaledWeights[!is.finite(scaledWeights)] = 0

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

makeZizMarginalPdf = function(posteriorGrid) {
  list(
    pi = approxfun(
      posteriorGrid$pi,
      posteriorGrid$marginalDensity$pi,
      yleft = 0,
      yright = 0,
      rule = 1
    ),
    shape = approxfun(
      posteriorGrid$shape,
      posteriorGrid$marginalDensity$shape,
      yleft = 0,
      yright = 0,
      rule = 1
    )
  )
}

fitZIDistBayesNumerical = function(x,
                                   nterms = 10,
                                   prior = makePrior(),
                                   shape1 = 1,
                                   shape2 = 1,
                                   nPiGrid = 101,
                                   nShapeGrid = 101,
                                   ...) {
  nvals = 1:nterms
  posteriorGrid = makeZizPosteriorGrid(
    x = x,
    prior = prior,
    shape1 = shape1,
    shape2 = shape2,
    nPiGrid = nPiGrid,
    nShapeGrid = nShapeGrid
  )

  par = posteriorGrid$mean
  marginalPdf = makeZizMarginalPdf(posteriorGrid)

  fitted = (1 - par[["pi"]]) * dzetaStandard(nvals, shape = par[["shape"]])
  fitted[nvals == 1] = fitted[nvals == 1] + par[["pi"]]
  names(fitted) = if (x$type == "P") {
    paste0("P", nvals - 1)
  } else {
    paste0("S", nvals)
  }

  result = list(
    psData = x,
    fit = list(par = par),
    pi = unname(par[["pi"]]),
    shape = unname(par[["shape"]]),
    var.cov = posteriorGrid$varCov,
    fitted = fitted,
    posteriorGrid = posteriorGrid,
    marginalPdf = marginalPdf,
    pdf = marginalPdf$shape,
    model = "ziz",
    method = "bayes",
    posteriorMethod = "numerical"
  )

  class(result) = "psFit"
  result
}
