test_that("Bayesian zero-inflated fits contain psPosterior objects", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(8, 3, 1),
    type = "P"
  )

  fit = fitZIDist(
    pData,
    nterms = 3,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "numerical"),
    nPiGrid = 21,
    nShapeGrid = 21
  )

  expect_s3_class(fit, "psFit")
  expect_s3_class(fit$posterior, "psPosterior")
  expect_identical(fit$posterior$method, "numerical")
  expect_equal(fit$posterior$probabilities, fit$posteriorProbs)
  expect_identical(fit$posterior$representation, fit$posteriorRepresentation)
  expect_named(
    fit$posterior,
    c(
      "method", "parameters", "probabilities", "representation",
      "level", "diagnostics", "model"
    )
  )
})

test_that("psPosterior print summary and fitted methods are coherent", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(8, 3, 1),
    type = "P"
  )

  fit = fitZIDist(
    pData,
    nterms = 3,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "numerical"),
    nPiGrid = 21,
    nShapeGrid = 21
  )

  expect_output(print(fit$posterior), "fitPS posterior approximation")
  posteriorSummary = summary(fit$posterior)
  expect_s3_class(posteriorSummary, "summary.psPosterior")
  expect_output(print(posteriorSummary), "Summary of fitPS posterior")
  expect_s3_class(summary(fit), "summary.psPosterior")

  posteriorFitted = fitted(fit$posterior)
  expect_equal(
    unname(posteriorFitted),
    fit$posterior$probabilities$estimate
  )
  expect_identical(
    names(posteriorFitted),
    fit$posterior$probabilities$term
  )
})

test_that("all posterior engines attach their native representation", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(8, 3, 1),
    type = "P"
  )

  mcmcFit = fitZIDist(
    pData,
    nterms = 2,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "mcmc"),
    nIter = 1000,
    nBurnIn = 200,
    seed = 17
  )
  expect_s3_class(mcmcFit$posterior, "psPosterior")
  expect_identical(mcmcFit$posterior$representation, mcmcFit$posteriorRepresentation)

  importanceFit = fitZIDist(
    pData,
    nterms = 2,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "importance"),
    nSamples = 200,
    seed = 18
  )
  expect_s3_class(importanceFit$posterior, "psPosterior")
  expect_identical(
    importanceFit$posterior$representation,
    importanceFit$posteriorRepresentation
  )
  expect_identical(
    importanceFit$posterior$diagnostics,
    importanceFit$posteriorDiagnostics
  )
})

test_that("posteriorProbs dispatches for psFit and psPosterior objects", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(8, 3, 1),
    type = "P"
  )

  fit = fitZIDist(
    pData,
    nterms = 4,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "numerical"),
    nPiGrid = 21,
    nShapeGrid = 21
  )

  expect_equal(posteriorProbs(fit), fit$posterior$probabilities)
  expect_equal(posteriorProbs(fit$posterior), fit$posterior$probabilities)
  expect_equal(
    posteriorProbs(fit, n = 2)$term,
    c("P0", "P1")
  )
  expect_equal(
    posteriorProbs(fit, n = c(0, 2))$term,
    c("P0", "P2")
  )
})

test_that("fitted distinguishes plug-in and posterior mean probabilities", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(8, 3, 1),
    type = "P"
  )

  fit = fitZIDist(
    pData,
    nterms = 3,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "numerical"),
    nPiGrid = 21,
    nShapeGrid = 21
  )

  expect_equal(fitted(fit), fit$fitted)
  expect_equal(fitted(fit, type = "plugIn"), fit$fitted)
  expect_equal(
    unname(fitted(fit, type = "posteriorMean")),
    fit$posterior$probabilities$estimate
  )
  expect_identical(
    names(fitted(fit, type = "posteriorMean")),
    fit$posterior$probabilities$term
  )
})

test_that("posteriorProbs rejects non-Bayesian fits", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(8, 3, 1),
    type = "P"
  )

  fit = fitZIDist(pData, nterms = 3)

  expect_error(
    posteriorProbs(fit),
    "only available for Bayesian"
  )
  expect_error(
    fitted(fit, type = "posteriorMean"),
    "require a Bayesian"
  )
})

test_that("plain zeta Bayesian fits use the common psPosterior contract", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(8, 3, 1),
    type = "P"
  )
  prior = makePrior(family = "uniform", range = c(1.1, 4))

  numericalFit = fitDist(
    pData,
    nterms = 4,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "numerical", prior = prior)
  )
  set.seed(2301)
  mcmcFit = fitDist(
    pData,
    nterms = 4,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "mcmc", prior = prior),
    nIter = 1000,
    nBurnIn = 100,
    silent = TRUE
  )

  for (fit in list(numericalFit, mcmcFit)) {
    expect_s3_class(fit$posterior, "psPosterior")
    expect_identical(fit$posterior$representation, fit$posteriorRepresentation)
    expect_identical(fit$posterior$model, "zeta")
    expect_identical(fit$posteriorProbs, fit$posterior$probabilities)
    expect_identical(fit$posterior$parameters$parameter, "shape")
    expect_identical(fit$posterior$probabilities$term, paste0("P", 0:3))
    expect_null(dim(fit$fitted))
  }
})

test_that("Bayesian fit finalisation preserves fitted vector contracts", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(8, 3, 1),
    type = "P"
  )

  zizFit = fitZIDist(
    pData,
    nterms = 4,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "numerical"),
    nPiGrid = 21,
    nShapeGrid = 21
  )

  expect_null(dim(zizFit$fitted))
  expect_identical(names(zizFit$fitted), paste0("P", 0:3))
  expect_identical(
    zizFit$posterior$representation,
    zizFit$posteriorRepresentation
  )
})
