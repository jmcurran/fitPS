test_that("default Bayesian transforms provide a validated identity contract", {
  model = psModel(
    model = "identity",
    parameterNames = c("alpha", "beta"),
    subclass = "identityBayesModel"
  )
  parameters = c(alpha = -0.5, beta = 2)

  expect_identical(modelToWorking(model, parameters), parameters)
  expect_identical(modelFromWorking(model, parameters), parameters)
  expect_identical(modelWorkingLogJacobian(model, parameters), 0)

  expect_error(
    modelToWorking(model, c(alpha = 1)),
    "matching modelParameterNames"
  )
})


test_that("external models can supply prior mathematics and Bayesian controls", {
  model = psModel(
    model = "external-bayes",
    parameterNames = "theta",
    subclass = "externalBayesContractModel",
    supportedEngines = "mcmc"
  )
  x = makePSData(n = 0:2, count = c(10, 4, 1), type = "P")
  prior = list(mean = 0, sd = 2)

  modelLogPrior.externalBayesContractModel = function(model, parameters, prior, ...) {
    dnorm(parameters[["theta"]], mean = prior$mean, sd = prior$sd, log = TRUE)
  }
  modelBayesControl.externalBayesContractModel = function(model, x, engine, prior, ...) {
    list(start = c(theta = prior$mean))
  }

  namespace = asNamespace("fitPS")
  registerS3method(
    "modelLogPrior",
    "externalBayesContractModel",
    modelLogPrior.externalBayesContractModel,
    envir = namespace
  )
  registerS3method(
    "modelBayesControl",
    "externalBayesContractModel",
    modelBayesControl.externalBayesContractModel,
    envir = namespace
  )

  expect_equal(
    modelLogPrior(model, c(theta = 1), prior),
    dnorm(1, mean = 0, sd = 2, log = TRUE)
  )
  expect_identical(
    modelBayesControl(model, x, engine = "mcmc", prior = prior)$start,
    c(theta = 0)
  )
  expect_identical(supportedPosteriorEngines(model), "mcmc")
})


test_that("built-in models expose prior mathematics through the public contract", {
  prior = makePrior(family = "custom", range = c(1.1, 4), logd = function(shape) {
    dunif(shape, min = 1.1, max = 4, log = TRUE)
  })

  expect_equal(
    modelLogPrior(zetaModel(), c(shape = 2), prior),
    prior$logd(2)
  )
  expect_equal(
    modelLogPrior(
      zizModel(),
      c(pi = 0.25, shape = 2),
      prior,
      shape1 = 2,
      shape2 = 3
    ),
    dbeta(0.25, 2, 3, log = TRUE) + prior$logd(2)
  )

  logPrior = makePrior(family = "custom", range = c(0.1, 0.9), logd = function(pi) {
    dunif(pi, min = 0.1, max = 0.9, log = TRUE)
  })
  expect_equal(
    modelLogPrior(logarithmicModel(), c(pi = 0.4), logPrior),
    logPrior$logd(0.4)
  )
})


test_that("built-in working transformations round trip and report Jacobians", {
  zeta = zetaModel()
  zetaParameters = c(shape = 2.5)
  zetaWorking = modelToWorking(zeta, zetaParameters)
  expect_equal(modelFromWorking(zeta, zetaWorking), zetaParameters)
  expect_equal(modelWorkingLogJacobian(zeta, zetaWorking), unname(zetaWorking))

  logarithmic = logarithmicModel()
  logParameters = c(pi = 0.35)
  logWorking = modelToWorking(logarithmic, logParameters)
  expect_equal(modelFromWorking(logarithmic, logWorking), logParameters)
  expect_equal(
    modelWorkingLogJacobian(logarithmic, logWorking),
    log(0.35) + log1p(-0.35)
  )

  ziz = zizModel()
  zizParameters = c(pi = 0.3, shape = 2.4)
  zizWorking = modelToWorking(ziz, zizParameters)
  expect_equal(modelFromWorking(ziz, zizWorking), zizParameters)
  expect_equal(
    modelWorkingLogJacobian(ziz, zizWorking),
    log(0.3) + log1p(-0.3) + log(1.4)
  )
})


test_that("built-in Bayesian controls provide named natural-scale starts", {
  x = makePSData(n = 0:2, count = c(10, 4, 1), type = "P")
  shapePrior = makePrior(family = "custom", range = c(1.1, 4), logd = function(shape) {
    dunif(shape, min = 1.1, max = 4, log = TRUE)
  })
  logPrior = makePrior(family = "custom", range = c(0.1, 0.9), logd = function(pi) {
    dunif(pi, min = 0.1, max = 0.9, log = TRUE)
  })

  expect_identical(
    modelBayesControl(zetaModel(), x, "mcmc", shapePrior)$start,
    c(shape = 2)
  )
  expect_identical(
    modelBayesControl(logarithmicModel(), x, "mcmc", logPrior)$start,
    c(pi = 0.5)
  )
  expect_identical(
    modelBayesControl(zizModel(), x, "mcmc", shapePrior)$start,
    c(pi = 0.5, shape = 2)
  )
})
