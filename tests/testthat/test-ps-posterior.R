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
  expect_identical(fit$posterior$representation, fit$posteriorGrid)
  expect_named(
    fit$posterior,
    c(
      "method", "parameters", "probabilities", "representation",
      "level", "diagnostics"
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
  expect_identical(mcmcFit$posterior$representation, mcmcFit$chain)

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
    importanceFit$importance
  )
  expect_identical(
    importanceFit$posterior$diagnostics,
    importanceFit$posteriorDiagnostics
  )
})
