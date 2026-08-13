legacyNumericalZetaFit = function(x, prior, nterms) {
  nvals = seq_len(nterms)
  observations = if (x$type == "P") {
    x$data$n + 1L
  } else {
    x$data$n
  }

  logLikelihood = function(shape) {
    sum(x$data$rn * dzetaStandard(observations, shape = shape, log = TRUE))
  }
  negLogLikelihoodTimesPrior = function(shape) {
    -(logLikelihood(shape) + prior$logd(shape))
  }

  optimum = optim(
    par = mean(prior$range),
    fn = negLogLikelihoodTimesPrior,
    method = "L-BFGS-B",
    lower = prior$range[1],
    upper = prior$range[2]
  )
  logShift = optimum$value

  unnormalisedDensity = Vectorize(function(shape) {
    exp(-negLogLikelihoodTimesPrior(shape) + logShift)
  })
  normalisingIntegral = integrate(
    unnormalisedDensity,
    lower = prior$range[1],
    upper = prior$range[2]
  )
  logNormalisingConstant = log(normalisingIntegral$value)

  density = Vectorize(function(shape) {
    exp(
      -negLogLikelihoodTimesPrior(shape) +
        logShift -
        logNormalisingConstant
    )
  })
  meanIntegral = integrate(
    function(shape) {
      shape * density(shape)
    },
    lower = prior$range[1],
    upper = prior$range[2]
  )
  secondMomentIntegral = integrate(
    function(shape) {
      shape^2 * density(shape)
    },
    lower = prior$range[1],
    upper = prior$range[2]
  )

  shape = meanIntegral$value
  variance = secondMomentIntegral$value - shape^2
  fitted = dzetaStandard(nvals, shape = shape)
  names(fitted) = if (x$type == "P") {
    paste0("P", nvals - 1L)
  } else {
    paste0("S", nvals)
  }

  list(
    shape = shape,
    var.shape = variance,
    fitted = fitted,
    pdf = density
  )
}

test_that("zeta model log likelihood matches the direct likelihood", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")
  shape = 2.2
  observations = pData$data$n + 1L
  expected = sum(
    pData$data$rn * dzetaStandard(observations, shape = shape, log = TRUE)
  )

  actual = modelLogLikelihood(
    zetaModel(),
    parameters = list(shape = shape),
    data = pData
  )

  expect_equal(actual, expected)
})

test_that("numerical engine fits and summarises the plain zeta posterior", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 4))
  engine = numericalPosteriorEngine()
  model = zetaModel()

  representation = fitPosterior(engine, model, pData, prior)
  summary = summarisePosterior(engine, model, representation)
  pointEstimate = posteriorPointEstimate(engine, model, representation)
  diagnostics = posteriorDiagnostics(engine, representation)

  expect_s3_class(representation, "numericalPosteriorRepresentation")
  expect_true(is.function(representation$value$density))
  expect_named(summary, c("parameter", "estimate", "sd"))
  expect_identical(summary$parameter, "shape")
  expect_equal(unname(pointEstimate), summary$estimate)
  expect_named(pointEstimate, "shape")
  expect_true(summary$sd > 0)
  expect_identical(diagnostics$model, "zeta")
  expect_equal(
    integrate(representation$value$density, 1.1, 4)$value,
    1,
    tolerance = 1e-6
  )
})

test_that("Stage 6.4 numerical zeta fit reproduces the legacy algorithm", {
  prior = makePrior(family = "uniform", range = c(1.1, 4))

  for (type in c("P", "S")) {
    data = if (type == "P") {
      makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = type)
    } else {
      makePSData(n = c(1, 2, 3), count = c(8, 3, 1), type = type)
    }

    expected = legacyNumericalZetaFit(data, prior, nterms = 4)
    actual = fitDistBayesIntegrate(data, prior = prior, nterms = 4)

    expect_equal(actual$shape, expected$shape, tolerance = 1e-10)
    expect_equal(actual$var.shape, expected$var.shape, tolerance = 1e-10)
    expect_equal(actual$fitted, expected$fitted, tolerance = 1e-10)
    grid = seq(1.2, 3.8, length.out = 7)
    expect_equal(actual$pdf(grid), expected$pdf(grid), tolerance = 1e-10)
    expect_s3_class(
      actual$posteriorRepresentation,
      "numericalPosteriorRepresentation"
    )
  }
})

test_that("fitDist numerical Bayesian path uses the migrated engine", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 4))

  fit = fitDist(
    pData,
    method = "bayes",
    bayesOptions = list(
      posteriorMethod = "numerical",
      prior = prior
    ),
    nterms = 4
  )

  expect_identical(fit$method, "bayes")
  expect_identical(fit$posteriorMethod, "numerical")
  expect_s3_class(
    fit$posteriorRepresentation,
    "numericalPosteriorRepresentation"
  )
  expect_true(is.function(fit$pdf))
})

test_that("numerical ZIZ is migrated by Stage 6.6", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 4))

  representation = fitPosterior(
    numericalPosteriorEngine(),
    zizModel(),
    pData,
    prior,
    nPiGrid = 9,
    nShapeGrid = 9
  )

  expect_s3_class(representation, "numericalPosteriorRepresentation")
})
