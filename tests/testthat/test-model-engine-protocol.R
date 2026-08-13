test_that("model descriptors declare parameters and supported engines", {
  zeta = zetaModel()
  ziz = zizModel()

  expect_s3_class(zeta, "zetaModel")
  expect_s3_class(zeta, "psModel")
  expect_identical(modelParameterNames(zeta), "shape")
  expect_identical(supportedPosteriorEngines(zeta), c("numerical", "mcmc"))

  expect_s3_class(ziz, "zizModel")
  expect_s3_class(ziz, "psModel")
  expect_identical(modelParameterNames(ziz), c("pi", "shape"))
  expect_identical(
    supportedPosteriorEngines(ziz),
    c("numerical", "mcmc", "laplace", "importance")
  )
})

test_that("model probability dispatch reproduces direct zeta-family helpers", {
  zetaExpected = zetaProbabilities(c(2, 2.5), n = 0:2, type = "P")
  zetaActual = modelProbabilities(
    zetaModel(),
    parameters = list(shape = c(2, 2.5)),
    n = 0:2,
    type = "P"
  )
  expect_equal(zetaActual, zetaExpected)

  zizExpected = zizProbabilities(
    pi = c(0.2, 0.3),
    shape = c(2, 2.5),
    n = 1:3,
    type = "S"
  )
  zizActual = modelProbabilities(
    zizModel(),
    parameters = list(pi = c(0.2, 0.3), shape = c(2, 2.5)),
    n = 1:3,
    type = "S"
  )
  expect_equal(zizActual, zizExpected)
})

test_that("shared observation mapping is consistent for zeta and ZIZ models", {
  pData = makePSData(n = c(0, 1, 2), count = c(5, 3, 1), type = "P")
  sData = makePSData(n = c(1, 2, 3), count = c(5, 3, 1), type = "S")

  expect_identical(modelObservationData(zetaModel(), pData), c(1L, 2L, 3L))
  expect_identical(modelObservationData(zizModel(), pData), c(1L, 2L, 3L))
  expect_identical(modelObservationData(zetaModel(), sData), c(1L, 2L, 3L))
  expect_identical(zizObservationData(pData), psObservationData(pData))
})

test_that("probability helper validation and naming are shared", {
  expect_identical(psProbabilityTermNames(0:2, "P"), c("P0", "P1", "P2"))
  expect_identical(psProbabilityTermNames(1:3, "S"), c("S1", "S2", "S3"))
  expect_identical(latentPsValues(0:2, "P"), c(1L, 2L, 3L))
  expect_identical(latentPsValues(1:3, "S"), c(1L, 2L, 3L))

  expect_error(normaliseProbabilityIndices(-1, "P"), "non-negative")
  expect_error(normaliseProbabilityIndices(0, "S"), "positive")
})

test_that("posterior-engine descriptors expose the sparse support matrix", {
  numerical = numericalPosteriorEngine()
  mcmc = mcmcPosteriorEngine()
  laplace = laplacePosteriorEngine()
  importance = importancePosteriorEngine()

  expect_s3_class(numerical, "psPosteriorEngine")
  expect_identical(posteriorEngineName(numerical), "numerical")
  expect_identical(posteriorEngineName(mcmc), "mcmc")
  expect_identical(posteriorEngineName(laplace), "laplace")
  expect_identical(posteriorEngineName(importance), "importance")

  expect_true(supportsPosteriorEngine(zetaModel(), numerical))
  expect_true(supportsPosteriorEngine(zetaModel(), mcmc))
  expect_false(supportsPosteriorEngine(zetaModel(), laplace))
  expect_true(supportsPosteriorEngine(zizModel(), laplace))
  expect_true(supportsPosteriorEngine(zizModel(), importance))
})

test_that("posterior representation contract preserves engine-specific payloads", {
  engine = mcmcPosteriorEngine()
  chain = data.frame(shape = c(2, 2.1, 2.2))
  representation = newPsPosteriorRepresentation(
    engine,
    value = chain,
    metadata = list(seed = 123L)
  )

  expect_s3_class(representation, "mcmcPosteriorRepresentation")
  expect_s3_class(representation, "psPosteriorRepresentation")
  expect_identical(representation$method, "mcmc")
  expect_identical(representation$value, chain)
  expect_identical(representation$metadata$seed, 123L)
  expect_invisible(validatePosteriorRepresentation(representation, engine))

  expect_error(
    validatePosteriorRepresentation(representation, numericalPosteriorEngine()),
    "does not match"
  )
})

test_that("unmigrated and unsupported engine paths fail explicitly", {
  pData = makePSData(n = c(0, 1, 2), count = c(5, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 4))

  expect_error(
    fitPosterior(laplacePosteriorEngine(), zetaModel(), pData, prior),
    "not supported"
  )
  expect_s3_class(
    fitPosterior(numericalPosteriorEngine(), zetaModel(), pData, prior),
    "numericalPosteriorRepresentation"
  )
  expect_error(
    modelLogLikelihood(zetaModel(), list(shape = 2), c(1L, 2L)),
    "data must be an object of class psData"
  )
})
