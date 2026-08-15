test_that("built-in MLE models share the model-comparison contract", {
  data(Psurveys)

  fits = list(
    zeta = fitDist(Psurveys$roux),
    ziz = fitZIDist(Psurveys$roux),
    logarithmic = fitlogDist(Psurveys$roux)
  )

  logLikelihoods = lapply(fits, logLik)
  parameterCounts = vapply(logLikelihoods, attr, integer(1L), which = "df")
  sampleSizes = vapply(logLikelihoods, attr, numeric(1L), which = "nobs")

  expect_identical(parameterCounts, c(zeta = 1L, ziz = 2L, logarithmic = 1L))
  expect_length(unique(sampleSizes), 1L)
  expect_true(all(vapply(logLikelihoods, function(value) {
    is.finite(as.numeric(value))
  }, logical(1L))))
  expect_true(all(vapply(fits, function(fit) {
    is.finite(deviance(fit)) && is.finite(AIC(fit)) && is.finite(BIC(fit))
  }, logical(1L))))
})


test_that("built-in models declare only meaningful posterior engines", {
  expect_identical(
    supportedPosteriorEngines(zetaModel()),
    c("numerical", "mcmc")
  )
  expect_identical(
    supportedPosteriorEngines(zizModel()),
    c("numerical", "mcmc", "laplace", "importance")
  )
  expect_identical(
    supportedPosteriorEngines(logarithmicModel()),
    c("numerical", "mcmc")
  )
})
