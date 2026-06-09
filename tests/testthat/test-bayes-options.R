test_that("normaliseBayesOptions supplies numerical posterior defaults", {
  options = normaliseBayesOptions()

  expect_identical(options$posteriorMethod, "numerical")
  expect_s3_class(options$prior, "psPrior")
})

test_that("normaliseBayesOptions accepts explicit posterior method and prior", {
  prior = makePrior(family = "uniform", range = c(1.1, 4))
  options = normaliseBayesOptions(
    bayesOptions = list(
      posteriorMethod = "mcmc",
      prior = prior
    )
  )

  expect_identical(options$posteriorMethod, "mcmc")
  expect_identical(options$prior, prior)
})

test_that("normaliseBayesOptions rejects ambiguous method and prior settings", {
  prior = makePrior(family = "uniform", range = c(1.1, 4))

  expect_error(
    normaliseBayesOptions(bayesOptions = list(method = "numerical")),
    "posteriorMethod"
  )
  expect_error(
    normaliseBayesOptions(
      bayesOptions = list(prior = prior),
      prior = prior
    ),
    "either as prior or bayesOptions"
  )
})

test_that("zero-inflated transformed parameter helpers round trip", {
  theta = c(pi = 0.25, shape = 2.5)
  working = zizThetaToWorking(theta)
  thetaAgain = zizWorkingToTheta(working)

  expect_named(working, c("eta", "tau"))
  expect_equal(unname(thetaAgain), unname(theta))
  expect_true(is.finite(zizWorkingLogJacobian(working)))
})

test_that("fitDist uses numerical posterior as canonical Bayesian default", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 4))
  fit = fitDist(
    pData,
    method = "bayes",
    bayesOptions = list(prior = prior)
  )

  expect_s3_class(fit, "psFit")
  expect_identical(fit$method, "bayes")
  expect_identical(fit$posteriorMethod, "numerical")
  expect_true(is.function(fit$pdf))
})

test_that("fitZIDist keeps MCMC available through posteriorMethod", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")
  fit = fitZIDist(
    pData,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "mcmc"),
    nIter = 1000,
    nBurnIn = 10
  )

  expect_s3_class(fit, "psFit")
  expect_identical(fit$method, "bayes")
  expect_identical(fit$posteriorMethod, "mcmc")
  expect_s3_class(fit$chain, "data.frame")
})
