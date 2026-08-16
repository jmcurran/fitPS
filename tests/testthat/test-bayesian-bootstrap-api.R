test_that("fit exposes Rubin Bayesian Bootstrap inference", {
  pData = makePSData(
    n = c(0, 1, 2, 3),
    count = c(40, 12, 4, 1),
    type = "P"
  )

  result = fit(
    pData,
    model = zetaModel(),
    method = "bayesianBootstrap",
    B = 16,
    seed = 1114,
    nterms = 4,
    level = 0.9
  )

  expect_s3_class(result, "psBayesianBootstrap")
  expect_identical(result$method, "rubin")
  expect_identical(result$model, "zeta")
  expect_identical(result$seed, 1114L)
  expect_equal(result$level, 0.9)
  expect_equal(nrow(result$replicates$parameters), 16)
  expect_equal(nrow(result$replicates$probabilities), 16)
  expect_equal(
    result$diagnostics$nSuccessful + result$diagnostics$nFailed,
    16
  )
})

test_that("Bayesian Bootstrap remains distinct from parametric Bayesian fitting", {
  pData = makePSData(
    n = c(0, 1, 2, 3),
    count = c(30, 10, 4, 1),
    type = "P"
  )

  result = fit(
    pData,
    model = logarithmicModel(),
    method = "bayesianBootstrap",
    B = 10,
    seed = 1115,
    nterms = 3
  )

  expect_s3_class(result, "psBayesianBootstrap")
  expect_false(inherits(result, "psPosterior"))
  expect_false(inherits(result, "psBootstrap"))
  expect_identical(result$method, "rubin")
})

test_that("Bayesian Bootstrap summary retains requested public information", {
  pData = makePSData(
    n = c(0, 1, 2, 3),
    count = c(35, 12, 5, 2),
    type = "P"
  )
  result = fit(
    pData,
    model = zizModel(),
    method = "bayesianBootstrap",
    B = 12,
    seed = 1116,
    nterms = 5
  )

  resultSummary = summary(result, nterms = 3)

  expect_s3_class(resultSummary, "summary.psBayesianBootstrap")
  expect_identical(resultSummary$method, "rubin")
  expect_identical(resultSummary$model, "ziz")
  expect_identical(resultSummary$seed, 1116L)
  expect_equal(nrow(resultSummary$probabilities), 3)
  expect_identical(resultSummary$parameters$parameter, c("pi", "shape"))
  expect_identical(resultSummary$diagnostics, result$diagnostics)
})

test_that("Bayesian Bootstrap print methods identify Rubin inference", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(20, 7, 2),
    type = "P"
  )
  result = fit(
    pData,
    model = zetaModel(),
    method = "bayesianBootstrap",
    B = 8,
    seed = 1117,
    nterms = 3
  )

  printed = capture.output(print(result))
  summaryPrinted = capture.output(print(summary(result, nterms = 2)))

  expect_true(any(grepl("Rubin Bayesian Bootstrap", printed, fixed = TRUE)))
  expect_true(any(grepl("Model: zeta", printed, fixed = TRUE)))
  expect_true(any(grepl("Rubin Bayesian Bootstrap", summaryPrinted, fixed = TRUE)))
  expect_true(any(grepl("Model-implied probability summaries", summaryPrinted, fixed = TRUE)))
})

test_that("fit validates Bayesian Bootstrap inputs", {
  pData = makePSData(n = c(0, 1, 2), count = c(10, 4, 1), type = "P")

  expect_error(
    fit(pData, model = "zeta", method = "bayesianBootstrap", B = 5),
    "model must inherit from psModel"
  )
  expect_error(
    fit(pData, model = zetaModel(), method = "bayesianBootstrap", B = 0),
    "B must be one positive integer"
  )
  expect_error(
    fit(
      pData,
      model = zetaModel(),
      method = "bayesianBootstrap",
      B = 5,
      level = 1
    ),
    "level must be one finite number strictly between 0 and 1"
  )

  result = fit(
    pData,
    model = zetaModel(),
    method = "bayesianBootstrap",
    B = 5,
    seed = 1118,
    nterms = 3
  )
  expect_error(summary(result, nterms = 0), "nterms must be one positive integer")
  expect_error(print(result, nterms = 0), "nterms must be one positive integer")
})


test_that("standalone bayesianBootstrap wrapper is no longer exported", {
  expect_false("bayesianBootstrap" %in% getNamespaceExports("fitPS"))
})
