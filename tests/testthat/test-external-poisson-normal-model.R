test_that("external Poisson-normal methods are registered for fitPS generics", {
  expect_identical(
    getS3method("modelObservationData", "externalPoissonNormalModel"),
    modelObservationData.externalPoissonNormalModel
  )
  expect_identical(
    getS3method("modelMleControl", "externalPoissonNormalModel"),
    modelMleControl.externalPoissonNormalModel
  )
  expect_identical(
    getS3method("modelProbabilities", "externalPoissonNormalModel"),
    modelProbabilities.externalPoissonNormalModel
  )
  expect_identical(
    getS3method("modelLogLikelihood", "externalPoissonNormalModel"),
    modelLogLikelihood.externalPoissonNormalModel
  )
})


test_that("external Poisson-normal probabilities integrate over a latent normal", {
  model = externalPoissonNormalModel()
  probabilities = modelProbabilities(
    model,
    parameters = list(mu = 0.2, sigma = 0.45),
    n = 0:30,
    type = "P"
  )

  expect_identical(modelParameterNames(model), c("mu", "sigma"))
  expect_equal(sum(probabilities), 1, tolerance = 1e-6)
  expect_true(all(probabilities >= 0))
})


test_that("external Poisson-normal preserves its probability shape across P and S surveys", {
  model = externalPoissonNormalModel()
  parameters = list(mu = 0.2, sigma = 0.45)

  pProbabilities = modelProbabilities(
    model,
    parameters = parameters,
    n = 0:8,
    type = "P"
  )
  sProbabilities = modelProbabilities(
    model,
    parameters = parameters,
    n = 1:9,
    type = "S"
  )

  expect_equal(as.numeric(pProbabilities), as.numeric(sProbabilities), tolerance = 1e-10)
})


test_that("external Poisson-normal maps matched P and S surveys to the same support", {
  counts = c(3032, 3240, 2035, 997, 426, 168, 64, 24, 9)
  pData = makePSData(n = 0:8, count = counts, type = "P")
  sData = makePSData(n = 1:9, count = counts, type = "S")
  model = externalPoissonNormalModel()
  parameters = list(mu = 0.2, sigma = 0.45)

  expect_equal(modelObservationData(model, pData), 0:8)
  expect_equal(modelObservationData(model, sData), 0:8)
  expect_equal(
    modelMleControl(model, pData)$start,
    modelMleControl(model, sData)$start
  )
  expect_equal(
    modelLogLikelihood(model, parameters, pData),
    modelLogLikelihood(model, parameters, sData),
    tolerance = 1e-8
  )
})


test_that("external Poisson-normal model fits through the generic MLE contract", {
  pData = makePSData(
    n = 0:8,
    count = c(3032, 3240, 2035, 997, 426, 168, 64, 24, 9),
    type = "P"
  )
  model = externalPoissonNormalModel()
  object = fit(pData, model = model, nterms = 9)

  expect_s3_class(object, "psFit")
  expect_s3_class(object$modelObject, "externalPoissonNormalModel")
  expect_identical(object$model, "poissonNormal")
  expect_true(is.finite(object$mu))
  expect_true(is.finite(object$sigma))
  expect_gt(object$sigma, 0)
  expect_equal(object$mu, 0.2, tolerance = 0.08)
  expect_equal(object$sigma, 0.45, tolerance = 0.08)
})


test_that("external Poisson-normal model participates in fitted and prediction APIs", {
  pData = makePSData(
    n = 0:8,
    count = c(3032, 3240, 2035, 997, 426, 168, 64, 24, 9),
    type = "P"
  )
  object = fit(pData, model = externalPoissonNormalModel(), nterms = 9)
  requested = 0:8
  expected = externalPoissonNormalProbability(
    requested,
    mu = object$mu,
    sigma = object$sigma
  )

  expect_equal(unname(fitted(object)), expected, tolerance = 1e-7)
  expect_equal(
    unname(predict(object, newdata = requested, interval = "none")),
    expected,
    tolerance = 1e-7
  )
})


test_that("external Poisson-normal model participates in MLE comparison criteria", {
  pData = makePSData(
    n = 0:8,
    count = c(3032, 3240, 2035, 997, 426, 168, 64, 24, 9),
    type = "P"
  )
  object = fit(pData, model = externalPoissonNormalModel(), nterms = 9)
  logLikelihood = logLik(object)

  expect_true(is.finite(as.numeric(logLikelihood)))
  expect_identical(attr(logLikelihood, "df"), 2L)
  expect_identical(attr(logLikelihood, "nobs"), sum(pData$data$rn))
  expect_equal(deviance(object), -2 * as.numeric(logLikelihood))
  expect_true(is.finite(AIC(object)))
  expect_true(is.finite(BIC(object)))
})


test_that("external Poisson-normal fits survive serialization", {
  pData = makePSData(
    n = 0:8,
    count = c(3032, 3240, 2035, 997, 426, 168, 64, 24, 9),
    type = "P"
  )
  object = fit(pData, model = externalPoissonNormalModel(), nterms = 9)
  path = tempfile(fileext = ".rds")
  saveRDS(object, path)
  restored = readRDS(path)

  expect_s3_class(restored$modelObject, "externalPoissonNormalModel")
  expect_equal(restored$mu, object$mu)
  expect_equal(restored$sigma, object$sigma)
  expect_equal(logLik(restored), logLik(object))
  expect_equal(fitted(restored), fitted(object))
})
