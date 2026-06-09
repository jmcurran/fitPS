test_that("default zeta prior is on standard shape scale", {
  prior = makePrior()

  expect_s3_class(prior, "psPrior")
  expect_gt(prior$range[1], 1)
  expect_true(is.finite(prior$logd(mean(prior$range))))
  expect_equal(prior$logd(1), -Inf)
})

test_that("makePrior rejects ranges that include invalid standard shape values", {
  expect_error(makePrior(range = c(1, 10)), "greater than 1")
  expect_error(makePrior(family = "uniform", range = c(0.5, 10)), "greater than 1")
})
