test_that("zero-inflated zeta MCMC is reproducible with a seed", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 4))

  fit1 = fitZIDist(
    pData,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "mcmc", prior = prior),
    shape1 = 2,
    shape2 = 5,
    theta0 = c(0.3, 2),
    nIter = 1200,
    nBurnIn = 300,
    seed = 779
  )
  fit2 = fitZIDist(
    pData,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "mcmc", prior = prior),
    shape1 = 2,
    shape2 = 5,
    theta0 = c(0.3, 2),
    nIter = 1200,
    nBurnIn = 300,
    seed = 779
  )

  expect_identical(fit1$chain, fit2$chain)
  expect_identical(fit1$fit$acceptance, fit2$fit$acceptance)
  expect_equal(fit1$pi, fit2$pi)
  expect_equal(fit1$shape, fit2$shape)
})

test_that("MCMC agrees with the numerical posterior under a non-uniform beta prior", {
  pData = makePSData(n = c(0, 1, 2, 3), count = c(12, 5, 2, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 4))

  numericalFit = fitZIDist(
    pData,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "numerical", prior = prior),
    shape1 = 2,
    shape2 = 5,
    nPiGrid = 81,
    nShapeGrid = 81
  )
  mcmcFit = fitZIDist(
    pData,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "mcmc", prior = prior),
    shape1 = 2,
    shape2 = 5,
    theta0 = c(0.3, 2),
    nIter = 12000,
    nBurnIn = 2000,
    seed = 20260806
  )

  expect_equal(mcmcFit$pi, numericalFit$pi, tolerance = 0.04)
  expect_equal(mcmcFit$shape, numericalFit$shape, tolerance = 0.12)
  piVarianceDifference = abs(
    unname(mcmcFit$var.cov["pi", "pi"]) -
      unname(numericalFit$var.cov["pi", "pi"])
  )
  expect_lte(piVarianceDifference, 0.001)
  expect_equal(
    unname(mcmcFit$var.cov["shape", "shape"]),
    unname(numericalFit$var.cov["shape", "shape"]),
    tolerance = 0.08
  )
  expect_true(all(mcmcFit$fit$acceptance > 0))
  expect_true(all(mcmcFit$fit$acceptance < 1))
})

test_that("legacy MCMC shape bounds define the prior actually used", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")

  fit = fitZIDist(
    pData,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "mcmc"),
    a = -1,
    b = 1,
    theta0 = c(0.4, 2),
    nIter = 1000,
    nBurnIn = 200,
    seed = 31
  )

  expect_equal(fit$bayesOptions$prior$range, 1 + exp(c(-1, 1)))
  expect_true(all(fit$chain$shape > fit$bayesOptions$prior$range[1]))
  expect_true(all(fit$chain$shape < fit$bayesOptions$prior$range[2]))
})
