test_that("MLE log likelihood uses the shared model contract", {
  data(Psurveys)
  fit = fitDist(Psurveys$roux)
  model = zetaModel()
  expected = modelLogLikelihood(model, c(shape = fit$shape), fit$psData)

  expect_equal(as.numeric(logLik(fit)), expected, tolerance = 1e-8)
  expect_equal(attr(logLik(fit), "df"), 1L)
  expect_equal(attr(logLik(fit), "nobs"), sum(fit$psData$data$rn))
  expect_equal(deviance(fit), -2 * expected, tolerance = 1e-8)
})

test_that("AIC and BIC distinguish zeta and ZIZ parameter counts", {
  data(Psurveys)
  zetaFit = fitDist(Psurveys$roux)
  zizFit = fitZIDist(Psurveys$roux)

  expect_equal(attr(logLik(zetaFit), "df"), 1L)
  expect_equal(attr(logLik(zizFit), "df"), 2L)
  expect_equal(AIC(zetaFit), -2 * as.numeric(logLik(zetaFit)) + 2)
  expect_equal(
    BIC(zizFit),
    -2 * as.numeric(logLik(zizFit)) + log(attr(logLik(zizFit), "nobs")) * 2
  )
})

test_that("AIC and BIC reject Bayesian fits", {
  data(Psurveys)
  fit = fitDist(Psurveys$roux, method = "bayes")

  expect_error(logLik(fit), "maximum-likelihood")
  expect_error(AIC(fit), "maximum-likelihood")
  expect_error(BIC(fit), "maximum-likelihood")
})

test_that("DIC uses numerical posterior deviance expectations", {
  data(Psurveys)
  fit = fitDist(Psurveys$roux, method = "bayes")
  value = DIC(fit)

  expect_true(is.finite(value))
  expect_true(is.finite(attr(value, "Dbar")))
  expect_true(is.finite(attr(value, "Dhat")))
  expect_equal(
    as.numeric(value),
    2 * attr(value, "Dbar") - attr(value, "Dhat")
  )
  expect_equal(
    attr(value, "pD"),
    attr(value, "Dbar") - attr(value, "Dhat")
  )
})

test_that("DIC rejects non-Bayesian and unsupported Laplace fits", {
  data(Psurveys)
  expect_error(DIC(fitDist(Psurveys$roux)), "Bayesian")

  fit = fitZIDist(
    Psurveys$roux,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "laplace")
  )
  expect_error(DIC(fit), "Laplace")
})
