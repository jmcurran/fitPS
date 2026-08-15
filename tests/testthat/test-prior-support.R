test_that("uniform priors can represent non-zeta parameter domains", {
  prior = makePrior(family = "uniform", range = c(0.1, 0.9))

  expect_s3_class(prior, "psPrior")
  expect_identical(prior$range, c(0.1, 0.9))
  expect_true(is.finite(prior$logd(0.5)))
  expect_identical(prior$logd(0.05), -Inf)
})

test_that("loguniform priors retain zeta shape-domain safeguards", {
  expect_error(
    makePrior(family = "loguniform", range = c(0.1, 0.9)),
    "greater than 1"
  )
})
