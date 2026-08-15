test_that("default zeta prior is on standard shape scale", {
  prior = makePrior()

  expect_s3_class(prior, "psPrior")
  expect_gt(prior$range[1], 1)
  expect_true(is.finite(prior$logd(mean(prior$range))))
  expect_equal(prior$logd(1), -Inf)
})

test_that("zeta posterior fitting rejects prior support below shape one", {
  pData = makePSData(n = c(0, 1, 2), count = c(5, 3, 1), type = "P")

  expect_error(makePrior(range = c(1, 10)), "greater than 1")
  expect_error(
    fitDist(
      pData,
      method = "bayes",
      prior = makePrior(family = "uniform", range = c(0.5, 10))
    ),
    "greater than 1"
  )
})
