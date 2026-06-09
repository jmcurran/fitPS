test_that("user-facing zeta shape must be greater than one", {
  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")

  expect_error(fitDist(pData, nterms = 4, start = 1), "greater than 1")
  expect_error(fitDist(pData, nterms = 4, shape = 1), "greater than 1")
  expect_error(rzeta(n = 1, shape = 1), "greater than 1")
  expect_error(rZIzeta(n = 1, pi = 0.5, shape = 1), "greater than 1")
})

test_that("probfun uses standard zeta shape for zeta P terms", {
  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")
  fit = fitDist(pData, nterms = 4, start = 2)
  fit$shape = 2

  observed = unname(probfun(fit)(0:2))
  expected = VGAM::dzeta(1:3, shape = 1)

  expect_equal(observed, expected)
})

test_that("probfun uses standard zeta shape for zeta S terms", {
  sData = makePSData(n = c(1, 2, 3), count = c(20, 5, 1), type = "S")
  fit = fitDist(sData, nterms = 4, start = 2)
  fit$shape = 2

  observed = unname(probfun(fit)(1:3))
  expected = VGAM::dzeta(1:3, shape = 1)

  expect_equal(observed, expected)
})

test_that("probfun uses standard zeta shape for zero-inflated P terms", {
  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")
  fit = fitZIDist(pData, nterms = 4, start = c(0.25, 2))
  fit$pi = 0.2
  fit$shape = 2

  observed = unname(probfun(fit)(0:2))
  expected = 0.8 * VGAM::dzeta(1:3, shape = 1)
  expected[1] = expected[1] + 0.2

  expect_equal(observed, expected)
})

test_that("predict.psFit uses standard zeta shape for new predictions", {
  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")
  fit = fitDist(pData, nterms = 4, start = 2)
  fit$shape = 2

  observed = unname(predict(fit, newdata = 0:2))
  expected = VGAM::dzeta(1:3, shape = 1)

  expect_equal(observed, expected)
})

test_that("predict.psFit uses standard zeta shape for new zero-inflated predictions", {
  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")
  fit = fitZIDist(pData, nterms = 4, start = c(0.25, 2))
  fit$pi = 0.2
  fit$shape = 2

  observed = unname(predict(fit, newdata = 0:2))
  expected = 0.8 * VGAM::dzeta(1:3, shape = 1)
  expected[1] = expected[1] + 0.2

  expect_equal(observed, expected)
})

test_that("fitted.psFit recomputes requested fitted values with standard shape", {
  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")
  fit = fitDist(pData, nterms = 4, start = 2)
  fit$shape = 2

  observed = unname(fitted(fit, n = 3))
  expected = VGAM::dzeta(1:3, shape = 1)

  expect_equal(observed, expected)
})

test_that("fitted psFit shape is stored as standard alpha", {
  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")
  fit = fitDist(pData, nterms = 4, start = 2)

  expect_gt(fit$shape, 1)
  expect_equal(unname(fit$fitted[1:3]), VGAM::dzeta(1:3, shape = fit$shape - 1))
})
