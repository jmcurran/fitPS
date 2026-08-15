test_that("public model API exports the external extension contract", {
  exports = getNamespaceExports("fitPS")
  expected = c(
    "psModel",
    "modelParameterNames",
    "modelObservationData",
    "modelProbabilities",
    "modelLogLikelihood",
    "modelMleControl",
    "supportedPosteriorEngines"
  )

  expect_true(all(expected %in% exports))
})


test_that("external Poisson S3 methods are registered for fitPS generics", {
  expect_identical(
    getS3method("modelObservationData", "externalPoissonModel"),
    modelObservationData.externalPoissonModel
  )
  expect_identical(
    getS3method("modelProbabilities", "externalPoissonModel"),
    modelProbabilities.externalPoissonModel
  )
  expect_identical(
    getS3method("modelLogLikelihood", "externalPoissonModel"),
    modelLogLikelihood.externalPoissonModel
  )
})


test_that("external Poisson preserves its probability shape across P and S surveys", {
  model = externalPoissonModel()
  parameters = list(lambda = 1.4)

  pProbabilities = modelProbabilities(
    model,
    parameters = parameters,
    n = 0:5,
    type = "P"
  )
  sProbabilities = modelProbabilities(
    model,
    parameters = parameters,
    n = 1:6,
    type = "S"
  )

  expect_equal(as.numeric(pProbabilities), as.numeric(sProbabilities))
  expect_equal(as.numeric(pProbabilities), dpois(0:5, lambda = 1.4))
})


test_that("external Poisson maps matched P and S surveys to the same support", {
  pData = makePSData(
    n = c(0, 1, 2, 3),
    count = c(30, 12, 5, 1),
    type = "P"
  )
  sData = makePSData(
    n = c(1, 2, 3, 4),
    count = c(30, 12, 5, 1),
    type = "S"
  )

  model = externalPoissonModel()
  expect_equal(modelObservationData(model, pData), 0:3)
  expect_equal(modelObservationData(model, sData), 0:3)

  pFit = fit(pData, model = model, nterms = 5)
  sFit = fit(sData, model = model, nterms = 5)
  expect_equal(pFit$lambda, sFit$lambda, tolerance = 1e-8)
  expect_equal(unname(fitted(pFit)), unname(fitted(sFit)), tolerance = 1e-8)
  expect_equal(as.numeric(logLik(pFit)), as.numeric(logLik(sFit)), tolerance = 1e-8)
})


test_that("external Poisson model fits through the public generic interface", {
  pData = makePSData(
    n = c(0, 1, 2, 3),
    count = c(30, 12, 5, 1),
    type = "P"
  )
  model = externalPoissonModel()
  object = fit(pData, model = model, nterms = 5)
  expectedLambda = weighted.mean(pData$data$n, pData$data$rn)

  expect_s3_class(object, "psFit")
  expect_s3_class(object$modelObject, "externalPoissonModel")
  expect_identical(object$model, "poisson")
  expect_equal(object$lambda, expectedLambda, tolerance = 1e-5)
  expect_identical(modelParameterNames(model), "lambda")
  expect_identical(supportedPosteriorEngines(model), c("numerical", "mcmc"))
})


test_that("external Poisson model participates in fitted and prediction APIs", {
  pData = makePSData(
    n = c(0, 1, 2, 3),
    count = c(30, 12, 5, 1),
    type = "P"
  )
  object = fit(pData, model = externalPoissonModel(), nterms = 5)
  requested = 0:4
  expected = dpois(requested, lambda = object$lambda)

  expect_equal(unname(fitted(object)), expected, tolerance = 1e-8)
  expect_equal(
    unname(predict(object, newdata = requested, interval = "none")),
    expected,
    tolerance = 1e-8
  )
})


test_that("external Poisson model participates in MLE comparison criteria", {
  pData = makePSData(
    n = c(0, 1, 2, 3),
    count = c(30, 12, 5, 1),
    type = "P"
  )
  object = fit(pData, model = externalPoissonModel(), nterms = 5)
  logLikelihood = logLik(object)

  expect_true(is.finite(as.numeric(logLikelihood)))
  expect_identical(attr(logLikelihood, "df"), 1L)
  expect_identical(attr(logLikelihood, "nobs"), sum(pData$data$rn))
  expect_equal(deviance(object), -2 * as.numeric(logLikelihood))
  expect_true(is.finite(AIC(object)))
  expect_true(is.finite(BIC(object)))
})


test_that("external model fits survive serialization when methods remain available", {
  pData = makePSData(
    n = c(0, 1, 2, 3),
    count = c(30, 12, 5, 1),
    type = "P"
  )
  object = fit(pData, model = externalPoissonModel(), nterms = 5)
  path = tempfile(fileext = ".rds")
  saveRDS(object, path)
  restored = readRDS(path)

  expect_s3_class(restored$modelObject, "externalPoissonModel")
  expect_equal(logLik(restored), logLik(object))
  expect_equal(fitted(restored), fitted(object))
})
