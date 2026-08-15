test_that("external Poisson fits through generic MCMC Bayesian orchestration", {
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
    nIter = 300L,
    nBurnIn = 100L,
    proposalScale = 0.2,
    seed = 779
  )

  expect_s3_class(fitObject, "psFit")
  expect_s3_class(fitObject$modelObject, "externalPoissonModel")
  expect_identical(fitObject$method, "bayes")
  expect_identical(fitObject$posteriorMethod, "mcmc")
  expect_true(is.finite(fitObject$lambda))
  expect_gt(fitObject$lambda, 0)
  expect_identical(fitObject$posterior$diagnostics$generic, TRUE)
  expect_equal(nrow(fitObject$posterior$representation$value$chain), 300L)
  expect_true(all(is.finite(fitted(fitObject))))
  expect_true(is.finite(DIC(fitObject)))
})


test_that("external generic MCMC requires an explicit model-specific prior", {
  data = makePSData(
    n = c(0, 1, 2),
    count = c(20, 8, 2),
    type = "P"
  )

  expect_error(
    fit(
      data,
      model = externalPoissonModel(),
      method = "bayes",
      nIter = 50L,
      nBurnIn = 20L
    ),
    "require an explicit model-specific prior"
  )
})


test_that("unsupported external posterior engines fail deliberately", {
  data = makePSData(
    n = c(0, 1, 2),
    count = c(20, 8, 2),
    type = "P"
  )

  expect_error(
    fit(
      data,
      model = externalPoissonModel(),
      method = "bayes",
      prior = list(shape = 2, rate = 1),
      bayesOptions = list(posteriorMethod = "laplace")
    ),
    "mcmc"
  )
})


test_that("external Poisson-normal uses the generic transformed MCMC path", {
  data = makePSData(
    n = 0:4,
    count = c(150, 120, 55, 20, 5),
    type = "P"
  )
  prior = list(muMean = 0, muSd = 2, sigmaScale = 1)

  fitObject = fit(
    data,
    model = externalPoissonNormalModel(),
    nterms = 5,
    method = "bayes",
    prior = prior,
    nIter = 40L,
    nBurnIn = 15L,
    proposalScale = c(mu = 0.08, sigma = 0.08),
    seed = 779
  )

  expect_s3_class(fitObject$modelObject, "externalPoissonNormalModel")
  expect_identical(fitObject$posteriorMethod, "mcmc")
  expect_true(is.finite(fitObject$mu))
  expect_true(is.finite(fitObject$sigma))
  expect_gt(fitObject$sigma, 0)

  chain = fitObject$posterior$representation$value$chain
  expect_identical(names(chain), c("mu", "sigma"))
  expect_true(all(chain$sigma > 0))
})
