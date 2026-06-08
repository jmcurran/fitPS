test_that("fitDist returns a zeta psFit for simple P data", {
  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")
  fit = fitDist(pData, nterms = 4)

  expect_s3_class(fit, "psFit")
  expect_identical(fit$model, "zeta")
  expect_identical(fit$method, "mle")
  expect_true(is.numeric(fit$shape))
  expect_length(fit$fitted, 4)
  expect_named(fit$fitted, paste0("P", 0:3))
})

test_that("predict.psFit returns named probabilities for fitted and new data", {
  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")
  fit = fitDist(pData, nterms = 4)

  fittedPredictions = predict(fit)
  newPredictions = predict(fit, newdata = 0:2)

  expect_named(fittedPredictions, paste0("P", 0:3))
  expect_named(newPredictions, paste0("P", 0:2))
  expect_true(all(is.finite(newPredictions)))
  expect_true(all(newPredictions >= 0))
  expect_true(all(newPredictions <= 1))
})

test_that("predict.psFit rejects non-integer new data", {
  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")
  fit = fitDist(pData, nterms = 4)

  expect_error(predict(fit, newdata = 1.5), "integers")
})

test_that("logLik.psFit returns a logLik object for fitted zeta models", {
  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")
  fit = fitDist(pData, nterms = 4)
  value = logLik(fit)

  expect_s3_class(value, "logLik")
  expect_true(is.finite(as.numeric(value)))
  expect_equal(attr(value, "df"), length(fit$fit$par))
  expect_equal(attr(value, "nobs"), sum(pData$data$rn))
})
