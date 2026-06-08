test_that("logLik.psFit returns a logLik object with model metadata", {
  psData = makePSData(n = 0:2, count = c(8, 3, 1), type = "P")
  fit = fitDist(psData, nterms = 3)

  ll = logLik(fit)

  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), -fit$fit$value)
  expect_equal(attr(ll, "df"), length(fit$fit$par))
  expect_equal(attr(ll, "nobs"), sum(psData$data$rn))
})

test_that("AIC and BIC dispatch through logLik.psFit", {
  psData = makePSData(n = 0:2, count = c(8, 3, 1), type = "P")
  fit = fitDist(psData, nterms = 3)
  ll = logLik(fit)

  expect_equal(AIC(fit), -2 * as.numeric(ll) + 2 * attr(ll, "df"))
  expect_equal(BIC(fit), -2 * as.numeric(ll) + log(attr(ll, "nobs")) * attr(ll, "df"))
})
