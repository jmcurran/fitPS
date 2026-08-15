test_that("ZIZ model log likelihood delegates to the characterised likelihood", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")
  model = zizModel()
  parameters = list(pi = 0.2, shape = 2.1)

  expected = zizLogLikelihood(
    obsData = zizObservationData(pData),
    counts = pData$data$rn,
    pi = parameters$pi,
    shape = parameters$shape
  )
  actual = modelLogLikelihood(model, parameters, pData)

  expect_equal(actual, expected, tolerance = 0)
})


test_that("numerical ZIZ fitting uses the generic two-dimensional cubature engine", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 4))
  engine = numericalPosteriorEngine()
  model = zizModel()

  representation = fitPosterior(
    engine,
    model,
    pData,
    prior,
    shape1 = 2,
    shape2 = 3,
    nPiGrid = 31,
    nShapeGrid = 31,
    tol = 1e-4
  )
  summary = summarisePosterior(engine, model, representation)
  pointEstimate = posteriorPointEstimate(engine, model, representation)
  diagnostics = posteriorDiagnostics(engine, representation)

  expect_s3_class(representation, "numericalPosteriorRepresentation")
  expect_identical(diagnostics$generic, TRUE)
  expect_identical(diagnostics$dimension, 2L)
  expect_identical(diagnostics$integrationMethod, "hcubature")
  expect_named(summary, c("parameter", "estimate", "sd"))
  expect_identical(summary$parameter, c("pi", "shape"))
  expect_equal(unname(pointEstimate), summary$estimate)
  expect_named(pointEstimate, c("pi", "shape"))
  expect_true(pointEstimate[["pi"]] > 0 && pointEstimate[["pi"]] < 1)
  expect_gt(pointEstimate[["shape"]], 1)
  expect_true(all(summary$sd >= 0))
  expect_true(is.finite(diagnostics$expectedDeviance))

  grid = representation$value$posteriorGrid
  expect_true(is.list(grid))
  expect_equal(grid$mean, pointEstimate, tolerance = 0)
  expect_equal(sum(grid$marginalDensity$pi) * grid$dPi, 1, tolerance = 0.08)
  expect_equal(sum(grid$marginalDensity$shape) * grid$dShape, 1, tolerance = 0.08)
})


test_that("numerical ZIZ fitting preserves matched P and S support mapping", {
  prior = makePrior(family = "uniform", range = c(1.1, 4))
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")
  sData = makePSData(n = c(1, 2, 3), count = c(8, 3, 1), type = "S")

  fitForData = function(data) {
    fit(
      data,
      model = zizModel(),
      nterms = 4,
      method = "bayes",
      prior = prior,
      bayesOptions = list(posteriorMethod = "numerical"),
      shape1 = 2,
      shape2 = 3,
      nPiGrid = 31,
      nShapeGrid = 31,
      tol = 1e-4
    )
  }

  pFit = fitForData(pData)
  sFit = fitForData(sData)

  expect_equal(c(pFit$pi, pFit$shape), c(sFit$pi, sFit$shape), tolerance = 1e-7)
  expect_equal(unname(pFit$fitted), unname(sFit$fitted), tolerance = 1e-7)
  expect_true(is.finite(DIC(pFit)))
  expect_true(is.finite(DIC(sFit)))
})
