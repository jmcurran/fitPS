test_that("fit accepts built-in model descriptors", {
  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")
  models = list(
    zeta = zetaModel(),
    ziz = zizModel(),
    logarithmic = logarithmicModel()
  )

  fits = lapply(models, function(model) {
    fit(pData, model = model, nterms = 4)
  })

  expect_true(all(vapply(fits, inherits, logical(1L), what = "psFit")))
  expect_identical(unname(vapply(fits, `[[`, character(1L), "model")), names(models))
  expect_true(all(vapply(fits, function(object) {
    inherits(object$modelObject, "psModel")
  }, logical(1L))))
  expect_true(all(vapply(fits, function(object) {
    is.finite(as.numeric(logLik(object)))
  }, logical(1L))))
})


test_that("established fitters retain their originating model object", {
  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")

  zetaFit = fitDist(pData, nterms = 4)
  zizFit = fitZIDist(pData, nterms = 4)
  logarithmicFit = fitlogDist(pData, nterms = 4)

  expect_s3_class(zetaFit$modelObject, "zetaModel")
  expect_s3_class(zizFit$modelObject, "zizModel")
  expect_s3_class(logarithmicFit$modelObject, "logarithmicModel")
})


test_that("model recovery prefers the retained model descriptor", {
  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")
  object = fitDist(pData, nterms = 4)
  retainedModel = object$modelObject
  object$model = "not-a-built-in-name"

  expect_identical(modelFromFit(object), retainedModel)
  expect_true(is.finite(as.numeric(logLik(object))))
})


test_that("legacy fits without modelObject retain built-in fallback behaviour", {
  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")
  object = fitDist(pData, nterms = 4)
  object$modelObject = NULL

  expect_s3_class(modelFromFit(object), "zetaModel")
  expect_true(is.finite(as.numeric(logLik(object))))
})


test_that("fit rejects objects that are not psModel descriptors", {
  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")

  expect_error(
    fit(pData, model = list(model = "zeta")),
    "model must inherit from psModel"
  )
})


test_that("fit forwards seeds to parametric Bayesian fitting", {
  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")
  options = list(posteriorMethod = "mcmc")

  fit1 = fit(
    pData,
    model = zizModel(),
    method = "bayes",
    bayesOptions = options,
    seed = 779,
    nterms = 4,
    nIter = 1020,
    nBurnIn = 20
  )
  fit2 = fit(
    pData,
    model = zizModel(),
    method = "bayes",
    bayesOptions = options,
    seed = 779,
    nterms = 4,
    nIter = 1020,
    nBurnIn = 20
  )

  expect_identical(
    fit1$posterior$representation$value$chain,
    fit2$posterior$representation$value$chain
  )
})
