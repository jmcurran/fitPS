test_that("Bayesian public fitters use the shared posterior contract", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(8, 3, 1),
    type = "P"
  )
  prior = makePrior(family = "uniform", range = c(1.1, 4))

  zetaFit = fitDist(
    pData,
    nterms = 3,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "numerical", prior = prior)
  )
  zizFit = fitZIDist(
    pData,
    nterms = 3,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "numerical", prior = prior),
    nPiGrid = 21,
    nShapeGrid = 21
  )

  for (fit in list(zetaFit, zizFit)) {
    expect_s3_class(fit, "psFit")
    expect_s3_class(fit$posterior, "psPosterior")
    expect_s3_class(fit$posterior$representation, "psPosteriorRepresentation")
    expect_false("posteriorRepresentation" %in% names(fit))
    expect_false("posteriorProbs" %in% names(fit))
  }
})

test_that("ZIZ Bayesian fits do not duplicate native posterior payloads at top level", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(8, 3, 1),
    type = "P"
  )

  numericalFit = fitZIDist(
    pData,
    nterms = 2,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "numerical"),
    nPiGrid = 21,
    nShapeGrid = 21
  )
  mcmcFit = fitZIDist(
    pData,
    nterms = 2,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "mcmc"),
    nIter = 1000,
    nBurnIn = 100,
    seed = 31
  )

  duplicatedFields = c(
    "posteriorGrid", "marginalPdf", "pdf", "chain", "laplace",
    "importance", "weightedSamples", "posteriorDiagnostics", "var.cov"
  )
  expect_false(any(duplicatedFields %in% names(numericalFit)))
  expect_false(any(duplicatedFields %in% names(mcmcFit)))
})

test_that("posteriorEngine constructs the requested strategy", {
  for (method in c("numerical", "mcmc", "laplace", "importance")) {
    engine = posteriorEngine(method)
    expect_s3_class(engine, "psPosteriorEngine")
    expect_identical(posteriorEngineName(engine), method)
  }
})
