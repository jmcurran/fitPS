test_that("fitDist honours the documented start argument", {
  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")

  fit = fitDist(pData, nterms = 4, start = 2)

  expect_s3_class(fit, "psFit")
  expect_identical(fit$model, "zeta")
  expect_true(is.finite(fit$shape))
})

test_that("fitDist preserves legacy shape start argument", {
  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")

  fit = fitDist(pData, nterms = 4, shape = 2)

  expect_s3_class(fit, "psFit")
  expect_identical(fit$model, "zeta")
  expect_true(is.finite(fit$shape))
})

test_that("fitZIDist honours the documented start argument", {
  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")

  fit = fitZIDist(pData, nterms = 4, start = c(0.25, 2))

  expect_s3_class(fit, "psFit")
  expect_identical(fit$model, "ziz")
  expect_true(is.finite(fit$pi))
  expect_true(is.finite(fit$shape))
})

test_that("fitZIDist preserves legacy shape start argument", {
  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")

  fit = fitZIDist(pData, nterms = 4, shape = c(0.25, 2))

  expect_s3_class(fit, "psFit")
  expect_identical(fit$model, "ziz")
  expect_true(is.finite(fit$pi))
  expect_true(is.finite(fit$shape))
})

test_that("predict.psFit uses fitted zero-inflation probability for ZIZ predictions", {
  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")
  fit = fitZIDist(pData, nterms = 4, start = c(0.25, 2))

  prediction = predict(fit, newdata = 0)
  expected = (1 - fit$pi) * VGAM::dzeta(1, shape = fit$shape - 1) + fit$pi

  expect_named(prediction, "P0")
  expect_equal(unname(prediction), expected)
  expect_true(prediction <= 1)
})
