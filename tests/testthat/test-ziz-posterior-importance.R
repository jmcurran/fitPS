test_that("zero-inflated importance sampler returns weighted posterior samples", {
  pData = makePSData(n = c(0, 1, 2), count = c(7, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 5))

  approximation = makeZizPosteriorImportance(
    x = pData,
    prior = prior,
    nSamples = 300,
    seed = 42,
    start = c(pi = 0.4, shape = 2)
  )

  expect_named(approximation$mean, c("pi", "shape"))
  expect_equal(nrow(approximation$samples), 300L)
  expect_equal(sum(approximation$samples$weight), 1, tolerance = 1e-10)
  expect_true(all(approximation$samples$weight >= 0))
  expect_true(approximation$mean[["pi"]] > 0)
  expect_true(approximation$mean[["pi"]] < 1)
  expect_true(approximation$mean[["shape"]] > 1)
  expect_true(is.finite(approximation$diagnostics$effectiveSampleSize))
  expect_true(approximation$diagnostics$effectiveSampleSize > 1)
  expect_true(approximation$diagnostics$maxWeight < 1)
  expect_equal(dim(approximation$varCov), c(2L, 2L))
})

test_that("zero-inflated importance sampler is reproducible with a seed", {
  pData = makePSData(n = c(0, 1, 2), count = c(7, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 5))

  approximation1 = makeZizPosteriorImportance(
    x = pData,
    prior = prior,
    nSamples = 300,
    seed = 123,
    start = c(pi = 0.4, shape = 2)
  )
  approximation2 = makeZizPosteriorImportance(
    x = pData,
    prior = prior,
    nSamples = 300,
    seed = 123,
    start = c(pi = 0.4, shape = 2)
  )

  expect_equal(approximation1$mean, approximation2$mean)
  expect_equal(approximation1$samples$weight, approximation2$samples$weight)
})

test_that("fitZIDist dispatches to Bayesian importance sampling", {
  pData = makePSData(n = c(0, 1, 2), count = c(7, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 5))

  fit = fitZIDist(
    pData,
    method = "bayes",
    bayesOptions = list(
      posteriorMethod = "importance",
      prior = prior
    ),
    nSamples = 300,
    seed = 99,
    start = c(pi = 0.4, shape = 2)
  )

  expect_s3_class(fit, "psFit")
  expect_equal(fit$method, "bayes")
  expect_equal(fit$posteriorMethod, "importance")
  expect_named(fit$fit$par, c("pi", "shape"))
  expect_equal(nrow(fit$weightedSamples), 300L)
  expect_true(is.list(fit$posteriorDiagnostics))
})
