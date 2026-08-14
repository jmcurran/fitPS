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
  expect_equal(nrow(fit$posterior$representation$value$approximation$samples), 300L)
  expect_true(is.list(fit$posterior$diagnostics))
})

test_that("importance posterior engine preserves ZIZ fit orchestration", {
  prior = makePrior(family = "uniform", range = c(1.1, 5))

  for (type in c("P", "S")) {
    x = if (type == "P") {
      makePSData(n = c(0, 1, 2), count = c(7, 3, 1), type = "P")
    } else {
      makePSData(n = c(1, 2, 3), count = c(7, 3, 1), type = "S")
    }

    approximation = makeZizPosteriorImportance(
      x = x,
      prior = prior,
      nSamples = 300,
      proposalScale = 2,
      seed = 271,
      start = c(pi = 0.4, shape = 2)
    )
    expectedPosteriorProbs = summariseZizSampleProbabilities(
      pi = approximation$samples$pi,
      shape = approximation$samples$shape,
      type = x$type,
      nterms = 4,
      weights = approximation$samples$weight,
      level = 0.95,
      posteriorMethod = "importance"
    )
    nvals = 1:4
    expectedFitted = (1 - approximation$mean[["pi"]]) *
      dzetaStandard(nvals, shape = approximation$mean[["shape"]])
    expectedFitted[nvals == 1] = expectedFitted[nvals == 1] +
      approximation$mean[["pi"]]
    names(expectedFitted) = if (type == "P") {
      paste0("P", nvals - 1)
    } else {
      paste0("S", nvals)
    }

    actual = fitZIDist(
      x = x,
      nterms = 4,
      prior = prior,
      method = "bayes",
      bayesOptions = list(posteriorMethod = "importance"),
      nSamples = 300,
      proposalScale = 2,
      seed = 271,
      start = c(pi = 0.4, shape = 2)
    )

    expect_s3_class(actual$posterior$representation, "importancePosteriorRepresentation")
    expect_equal(actual$fit$par, approximation$mean)
    expect_equal(actual$fitted, expectedFitted)
    expect_null(dim(actual$fitted))
    expect_equal(actual$posterior$probabilities, expectedPosteriorProbs)
  }
})

test_that("importance engine reports weighted-sample diagnostics", {
  pData = makePSData(n = c(0, 1, 2), count = c(7, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 5))

  representation = fitPosterior(
    importancePosteriorEngine(),
    zizModel(),
    pData,
    prior,
    nSamples = 300,
    seed = 19,
    start = c(pi = 0.4, shape = 2)
  )
  diagnostics = posteriorDiagnostics(importancePosteriorEngine(), representation)

  expect_s3_class(representation, "importancePosteriorRepresentation")
  expect_named(diagnostics, c("effectiveSampleSize", "maxWeight", "nSamples", "proposalScale"))
  expect_equal(diagnostics$nSamples, 300L)
  expect_error(
    fitPosterior(importancePosteriorEngine(), zetaModel(), pData, prior),
    "not supported"
  )
})
