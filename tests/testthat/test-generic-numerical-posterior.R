test_that("external one-parameter models fit through the generic numerical engine", {
  data = makePSData(
    n = c(0, 1, 2, 3),
    count = c(30, 12, 5, 1),
    type = "P"
  )
  prior = list(shape = 2, rate = 1)

  fitObject = fit(
    data,
    model = externalPoissonModel(),
    nterms = 5,
    method = "bayes",
    prior = prior,
    bayesOptions = list(posteriorMethod = "numerical")
  )

  posteriorShape = prior$shape + sum(data$data$rn * data$data$n)
  posteriorRate = prior$rate + sum(data$data$rn)
  expectedMean = posteriorShape / posteriorRate
  expectedVariance = posteriorShape / posteriorRate^2

  expect_s3_class(fitObject, "psFit")
  expect_s3_class(fitObject$modelObject, "externalPoissonModel")
  expect_identical(fitObject$posteriorMethod, "numerical")
  expect_equal(fitObject$lambda, expectedMean, tolerance = 1e-7)
  expect_equal(
    fitObject$posterior$representation$value$variance[[1L]],
    expectedVariance,
    tolerance = 1e-7
  )
  expect_identical(fitObject$posterior$diagnostics$generic, TRUE)
  expect_equal(
    integrate(
      fitObject$posterior$representation$value$density,
      lower = 0,
      upper = Inf
    )$value,
    1,
    tolerance = 1e-7
  )
  expect_true(all(is.finite(fitted(fitObject))))
  expect_true(is.finite(DIC(fitObject)))
})


test_that("built-in one-parameter numerical models use the generic engine", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")
  zetaPrior = makePrior(family = "uniform", range = c(1.1, 4))
  logarithmicPrior = makePrior(family = "uniform", range = c(0.1, 0.9))

  zetaRepresentation = fitPosterior(
    numericalPosteriorEngine(),
    zetaModel(),
    pData,
    zetaPrior
  )
  logarithmicRepresentation = fitPosterior(
    numericalPosteriorEngine(),
    logarithmicModel(),
    pData,
    logarithmicPrior
  )

  expect_identical(zetaRepresentation$metadata$generic, TRUE)
  expect_identical(logarithmicRepresentation$metadata$generic, TRUE)
  expect_equal(zetaRepresentation$metadata$bounds, c(lower = 1.1, upper = 4))
  expect_equal(
    logarithmicRepresentation$metadata$bounds,
    c(lower = 0.1, upper = 0.9)
  )
})


test_that("generic numerical fitting deliberately remains one-dimensional", {
  model = externalPoissonNormalModel()
  model$supportedEngines = c("numerical", model$supportedEngines)
  data = makePSData(n = 0:3, count = c(30, 15, 6, 2), type = "P")
  prior = list(muMean = 0, muSd = 2, sigmaScale = 1)

  expect_error(
    fit(
      data,
      model = model,
      method = "bayes",
      prior = prior,
      bayesOptions = list(posteriorMethod = "numerical")
    ),
    "one-parameter model"
  )
})
