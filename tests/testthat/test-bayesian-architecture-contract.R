test_that("Bayesian model and engine support matrix is explicit", {
  expect_identical(
    supportedPosteriorEngines(zetaModel()),
    c("numerical", "mcmc")
  )
  expect_identical(
    supportedPosteriorEngines(zizModel()),
    c("numerical", "mcmc", "laplace", "importance")
  )

  for (method in c("numerical", "mcmc")) {
    expect_true(supportsPosteriorEngine(zetaModel(), posteriorEngine(method)))
  }
  for (method in c("laplace", "importance")) {
    expect_false(supportsPosteriorEngine(zetaModel(), posteriorEngine(method)))
  }
  for (method in c("numerical", "mcmc", "laplace", "importance")) {
    expect_true(supportsPosteriorEngine(zizModel(), posteriorEngine(method)))
  }
})

test_that("provisional Bayesian fitters are absent from the package namespace", {
  namespace = asNamespace("fitPS")
  obsoleteFitters = c(
    "fitDistBayes",
    "fitDistBayesIntegrate",
    "fitZIDistBayes",
    "fitdistLaplace"
  )

  for (name in obsoleteFitters) {
    expect_false(exists(name, envir = namespace, inherits = FALSE))
  }
})

test_that("model descriptors can declare future distributions independently of engines", {
  alternativeModel = newPsModel(
    model = "alternative",
    parameterNames = "theta",
    supportedEngines = c("numerical", "mcmc"),
    subclass = "alternativeModel"
  )

  expect_s3_class(alternativeModel, "psModel")
  expect_identical(modelParameterNames(alternativeModel), "theta")
  expect_true(supportsPosteriorEngine(alternativeModel, numericalPosteriorEngine()))
  expect_true(supportsPosteriorEngine(alternativeModel, mcmcPosteriorEngine()))
  expect_false(supportsPosteriorEngine(alternativeModel, laplacePosteriorEngine()))
})
