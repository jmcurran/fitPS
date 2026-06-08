test_that("predict.psFit uses fitted zero-inflation probability for ZIZ P predictions", {
  psData = makePSData(n = 0:2, count = c(8, 3, 1), type = "P")
  fit = list(
    psData = psData,
    pi = 0.25,
    shape = 1,
    fitted = c(P0 = 0.625, P1 = 0.1875, P2 = 0.08333333),
    model = "ziz",
    method = "mle"
  )
  class(fit) = "psFit"

  prediction = predict(fit, newdata = 0)
  expected = (1 - fit$pi) * VGAM::dzeta(1, shape = fit$shape) + fit$pi

  expect_equal(unname(prediction), expected)
})

test_that("fitDist honours the start argument", {
  psData = makePSData(n = 0:2, count = c(8, 3, 1), type = "P")

  fit = fitDist(psData, nterms = 3, start = 0.5)

  expect_s3_class(fit, "psFit")
  expect_equal(fit$fit$par, fit$shape)
})

test_that("fitZIDist honours the start argument", {
  psData = makePSData(n = 0:2, count = c(8, 3, 1), type = "P")

  fit = fitZIDist(psData, nterms = 3, start = c(0.25, 1))

  expect_s3_class(fit, "psFit")
  expect_equal(fit$fit$par, c(fit$pi, fit$shape))
})
