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
    prior = prior
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
  expect_identical(
    fitObject$posterior$representation$value$grid$integrationRule,
    "simpson"
  )
  expect_equal(
    tail(fitObject$posterior$representation$value$grid$cumulative, 1L),
    1,
    tolerance = 1e-12
  )
  expect_true(all(diff(
    fitObject$posterior$representation$value$grid$cumulative
  ) >= 0))
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


test_that("external two-parameter models fit through adaptive cubature", {
  data = makePSData(n = 0:2, count = c(30, 12, 4), type = "P")
  prior = list(muMean = 0, muSd = 2, sigmaScale = 1)

  representation = fitPosterior(
    numericalPosteriorEngine(),
    externalPoissonNormalModel(),
    data,
    prior,
    tol = 1e-2,
    maxEval = 5000L,
    summaryGridSize = 11L
  )

  expect_s3_class(representation, "numericalPosteriorRepresentation")
  expect_identical(representation$metadata$dimension, 2L)
  expect_identical(representation$metadata$integrationMethod, "hcubature")
  expect_identical(representation$metadata$generic, TRUE)
  expect_named(representation$value$mean, c("mu", "sigma"))
  expect_true(all(is.finite(representation$value$mean)))
  expect_gt(representation$value$mean[["sigma"]], 0)
  expect_true(all(is.finite(representation$value$variance)))
  expect_true(is.finite(representation$metadata$expectedDeviance))
  expect_true(all(c("density", "mass", "integrationRule") %in% names(representation$value$grid)))
  expect_equal(sum(representation$value$grid$mass), 1, tolerance = 1e-12)
  expect_match(representation$value$grid$integrationRule, "simpson")
})


test_that("numerical fitting rejects models above two dimensions", {
  model = psModel(
    model = "three-parameter",
    parameterNames = c("alpha", "beta", "gamma"),
    supportedEngines = c("numerical", "mcmc"),
    subclass = "threeParameterNumericalModel"
  )
  data = makePSData(n = 0:2, count = c(20, 8, 2), type = "P")

  expect_error(
    fitPosterior(
      numericalPosteriorEngine(),
      model,
      data,
      prior = list()
    ),
    "at most two parameters.*MCMC"
  )
})


test_that("automatic Bayesian engine selection uses numerical only up to two dimensions", {
  expect_identical(defaultExternalPosteriorMethod(externalPoissonModel()), "numerical")
  expect_identical(defaultExternalPosteriorMethod(externalPoissonNormalModel()), "numerical")

  model = psModel(
    model = "three-parameter",
    parameterNames = c("alpha", "beta", "gamma"),
    supportedEngines = c("numerical", "mcmc"),
    subclass = "threeParameterAutomaticModel"
  )
  expect_identical(defaultExternalPosteriorMethod(model), "mcmc")
})
