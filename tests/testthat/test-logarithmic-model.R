test_that("logarithmic model declares its extension contract", {
  model = logarithmicModel()

  expect_s3_class(model, "logarithmicModel")
  expect_s3_class(model, "psModel")
  expect_identical(modelParameterNames(model), "pi")
  expect_identical(supportedPosteriorEngines(model), c("numerical", "mcmc"))
  expect_true(supportsPosteriorEngine(model, numericalPosteriorEngine()))
  expect_true(supportsPosteriorEngine(model, mcmcPosteriorEngine()))
  expect_false(supportsPosteriorEngine(model, laplacePosteriorEngine()))
  expect_false(supportsPosteriorEngine(model, importancePosteriorEngine()))
})

test_that("logarithmic probabilities use shared P and S support mapping", {
  pi = c(0.25, 0.6)
  pActual = modelProbabilities(
    logarithmicModel(),
    parameters = list(pi = pi),
    n = 0:2,
    type = "P"
  )
  sActual = modelProbabilities(
    logarithmicModel(),
    parameters = list(pi = pi),
    n = 1:3,
    type = "S"
  )

  pExpected = outer(pi, 1:3, Vectorize(function(parameter, value) {
    VGAM::dlog(value, shape = parameter)
  }))
  sExpected = pExpected
  colnames(pExpected) = c("P0", "P1", "P2")
  colnames(sExpected) = c("S1", "S2", "S3")

  expect_equal(pActual, pExpected)
  expect_equal(sActual, sExpected)
})

test_that("logarithmic likelihood follows the model contract for P and S data", {
  pData = makePSData(n = c(0, 1, 2), count = c(5, 3, 1), type = "P")
  sData = makePSData(n = c(1, 2, 3), count = c(5, 3, 1), type = "S")
  pi = 0.55

  pExpected = sum(pData$data$rn * VGAM::dlog(1:3, shape = pi, log = TRUE))
  sExpected = sum(sData$data$rn * VGAM::dlog(1:3, shape = pi, log = TRUE))

  expect_equal(
    modelLogLikelihood(logarithmicModel(), list(pi = pi), pData),
    pExpected
  )
  expect_equal(
    modelLogLikelihood(logarithmicModel(), list(pi = pi), sData),
    sExpected
  )
})

test_that("logarithmic MLE participates in common fitted-value and comparison APIs", {
  data(Psurveys)
  fit = fitlogDist(Psurveys$roux)

  expect_s3_class(fit, "psFit")
  expect_identical(fit$model, "logarithmic")
  expect_identical(fit$method, "mle")
  expect_true(fit$pi > 0 && fit$pi < 1)
  expect_true(is.finite(fit$var.pi) && fit$var.pi > 0)
  expect_equal(fitted(fit, n = 4), probfun(fit)(0:3))
  expect_equal(predict(fit, newdata = 0:3), probfun(fit)(0:3))
  expect_equal(attr(logLik(fit), "df"), 1L)
  expect_true(is.finite(AIC(fit)))
  expect_true(is.finite(BIC(fit)))
})

test_that("logarithmic MLE works for S data", {
  sData = makePSData(
    n = c(1, 2, 3, 4),
    count = c(20, 8, 3, 1),
    type = "S"
  )
  fit = fitlogDist(sData, nterms = 4)

  expect_identical(names(fit$fitted), c("S1", "S2", "S3", "S4"))
  expect_equal(fit$fitted, probfun(fit)(1:4))
})

test_that("logarithmic numerical posterior reuses the common posterior contract", {
  pData = makePSData(
    n = c(0, 1, 2, 3),
    count = c(20, 8, 3, 1),
    type = "P"
  )
  fit = fitlogDist(
    pData,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "numerical")
  )

  expect_s3_class(fit$posterior, "psPosterior")
  expect_identical(fit$posteriorMethod, "numerical")
  expect_identical(fit$posterior$parameters$parameter, "pi")
  expect_true(fit$pi > 0 && fit$pi < 1)
  expect_true(all(is.finite(fit$posterior$probabilities$estimate)))
  expect_true(is.finite(DIC(fit)))
})

test_that("logarithmic MCMC posterior supports summaries and fitted probabilities", {
  pData = makePSData(
    n = c(0, 1, 2, 3),
    count = c(20, 8, 3, 1),
    type = "P"
  )
  set.seed(728)
  fit = suppressWarnings(fitlogDist(
    pData,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "mcmc"),
    nIter = 200,
    nBurnIn = 100,
    silent = TRUE
  ))

  expect_identical(fit$posteriorMethod, "mcmc")
  expect_length(fit$posterior$representation$value$chain, 200L)
  expect_true(all(is.finite(fitted(fit, type = "posteriorMean"))))
  expect_true(is.finite(DIC(fit)))
})

test_that("logarithmic model rejects unsupported engines and invalid priors", {
  pData = makePSData(n = c(0, 1, 2), count = c(5, 3, 1), type = "P")

  expect_error(
    fitlogDist(
      pData,
      method = "bayes",
      bayesOptions = list(posteriorMethod = "laplace")
    ),
    "not supported"
  )
  expect_error(
    fitlogDist(
      pData,
      method = "bayes",
      prior = makePrior(family = "uniform", range = c(0.2, 1.2))
    ),
    "strictly inside"
  )
})
