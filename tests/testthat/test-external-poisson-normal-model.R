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


test_that("external Poisson-normal probabilities use mu and sigma", {
  model = externalPoissonNormalModel()
  probabilities = modelProbabilities(
    model,
    parameters = list(mu = 2, sigma = sqrt(3)),
    n = 0:20,
    type = "P"
  )

  expect_identical(modelParameterNames(model), c("mu", "sigma"))
  expect_equal(sum(probabilities), 1, tolerance = 1e-8)
  expect_true(all(probabilities >= 0))
})


test_that("external Poisson-normal model fits through the generic MLE contract", {
  pData = makePSData(
    n = 0:9,
    count = c(223, 223, 223, 149, 93, 48, 24, 10, 4, 2),
    type = "P"
  )
  model = externalPoissonNormalModel()
  object = fit(pData, model = model, nterms = 10)

  expect_s3_class(object, "psFit")
  expect_s3_class(object$modelObject, "externalPoissonNormalModel")
  expect_identical(object$model, "poissonNormal")
  expect_true(is.finite(object$mu))
  expect_true(is.finite(object$sigma))
  expect_gt(object$mu, 0)
  expect_gt(object$sigma, 0)
  expect_gte(object$sigma^2, object$mu)
  expect_lte(object$sigma^2, 2 * object$mu)
  expect_equal(
    object$mu,
    weighted.mean(pData$data$n, pData$data$rn),
    tolerance = 1e-3
  )
})


test_that("external Poisson-normal model participates in fitted and prediction APIs", {
  pData = makePSData(
    n = 0:9,
    count = c(223, 223, 223, 149, 93, 48, 24, 10, 4, 2),
    type = "P"
  )
  object = fit(pData, model = externalPoissonNormalModel(), nterms = 10)
  requested = 0:9
  expected = externalPoissonNormalProbability(
    requested,
    mu = object$mu,
    sigma = object$sigma
  )

  expect_equal(unname(fitted(object)), expected, tolerance = 1e-8)
  expect_equal(
    unname(predict(object, newdata = requested, interval = "none")),
    expected,
    tolerance = 1e-8
  )
})


test_that("external Poisson-normal model participates in MLE comparison criteria", {
  pData = makePSData(
    n = 0:9,
    count = c(223, 223, 223, 149, 93, 48, 24, 10, 4, 2),
    type = "P"
  )
  object = fit(pData, model = externalPoissonNormalModel(), nterms = 10)
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
    n = 0:9,
    count = c(223, 223, 223, 149, 93, 48, 24, 10, 4, 2),
    type = "P"
  )
  object = fit(pData, model = externalPoissonNormalModel(), nterms = 10)
  path = tempfile(fileext = ".rds")
  saveRDS(object, path)
  restored = readRDS(path)

  expect_s3_class(restored$modelObject, "externalPoissonNormalModel")
  expect_equal(restored$mu, object$mu)
  expect_equal(restored$sigma, object$sigma)
  expect_equal(logLik(restored), logLik(object))
  expect_equal(fitted(restored), fitted(object))
})
