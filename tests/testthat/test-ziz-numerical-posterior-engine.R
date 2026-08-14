legacyNumericalZizFit = function(x,
                                  prior,
                                  nterms,
                                  shape1 = 1,
                                  shape2 = 1,
                                  nPiGrid = 31,
                                  nShapeGrid = 31,
                                  level = 0.95) {
  nvals = seq_len(nterms)
  posteriorGrid = makeZizPosteriorGrid(
    x = x,
    prior = prior,
    shape1 = shape1,
    shape2 = shape2,
    nPiGrid = nPiGrid,
    nShapeGrid = nShapeGrid
  )

  par = posteriorGrid$mean
  posteriorProbs = summariseZizGridProbabilities(
    posteriorGrid = posteriorGrid,
    type = x$type,
    nterms = nterms,
    level = level
  )

  fitted = (1 - par[["pi"]]) * dzetaStandard(nvals, shape = par[["shape"]])
  fitted[nvals == 1] = fitted[nvals == 1] + par[["pi"]]
  names(fitted) = if (x$type == "P") {
    paste0("P", nvals - 1)
  } else {
    paste0("S", nvals)
  }

  list(
    par = par,
    var.cov = posteriorGrid$varCov,
    fitted = fitted,
    posteriorGrid = posteriorGrid,
    posteriorProbs = posteriorProbs
  )
}

test_that("ZIZ model log likelihood delegates to the characterised likelihood", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")
  model = zizModel()
  parameters = list(pi = 0.2, shape = 2.1)

  expected = zizLogLikelihood(
    obsData = zizObservationData(pData),
    counts = pData$data$rn,
    pi = parameters$pi,
    shape = parameters$shape
  )
  actual = modelLogLikelihood(model, parameters, pData)

  expect_equal(actual, expected, tolerance = 0)
})

test_that("numerical ZIZ engine implements the posterior protocol", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 4))
  engine = numericalPosteriorEngine()
  model = zizModel()

  representation = fitPosterior(
    engine,
    model,
    pData,
    prior,
    shape1 = 2,
    shape2 = 3,
    nPiGrid = 31,
    nShapeGrid = 31
  )
  summary = summarisePosterior(engine, model, representation)
  pointEstimate = posteriorPointEstimate(engine, model, representation)
  diagnostics = posteriorDiagnostics(engine, representation)

  expect_s3_class(representation, "numericalPosteriorRepresentation")
  expect_named(summary, c("parameter", "estimate", "sd"))
  expect_identical(summary$parameter, c("pi", "shape"))
  expect_equal(unname(pointEstimate), summary$estimate)
  expect_named(pointEstimate, c("pi", "shape"))
  expect_true(all(summary$sd >= 0))
  expect_identical(diagnostics$model, "ziz")
  expect_identical(diagnostics$nPiGrid, 31L)
  expect_identical(diagnostics$nShapeGrid, 31L)
  expect_equal(
    representation$value$posteriorGrid$mean,
    pointEstimate,
    tolerance = 0
  )
})

test_that("numerical ZIZ fitting reproduces the legacy orchestration", {
  prior = makePrior(family = "uniform", range = c(1.1, 4))

  for (type in c("P", "S")) {
    data = if (type == "P") {
      makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = type)
    } else {
      makePSData(n = c(1, 2, 3), count = c(8, 3, 1), type = type)
    }

    expected = legacyNumericalZizFit(
      data,
      prior,
      nterms = 4,
      shape1 = 2,
      shape2 = 3,
      nPiGrid = 31,
      nShapeGrid = 31
    )
    actual = fitZIDist(
      data,
      prior = prior,
      nterms = 4,
      method = "bayes",
      bayesOptions = list(posteriorMethod = "numerical"),
      shape1 = 2,
      shape2 = 3,
      nPiGrid = 31,
      nShapeGrid = 31
    )

    expect_equal(actual$fit$par, expected$par, tolerance = 0)
    expect_equal(actual$pi, unname(expected$par[["pi"]]), tolerance = 0)
    expect_equal(actual$shape, unname(expected$par[["shape"]]), tolerance = 0)
    expect_equal(actual$fitted, expected$fitted, tolerance = 1e-14)
    expect_equal(actual$posterior$probabilities, expected$posteriorProbs, tolerance = 0)
    expect_s3_class(
      actual$posterior$representation,
      "numericalPosteriorRepresentation"
    )
    expect_identical(actual$posterior$diagnostics$model, "ziz")
  }
})

test_that("fitZIDist numerical Bayesian path uses the migrated engine", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 4))

  fit = fitZIDist(
    pData,
    method = "bayes",
    bayesOptions = list(
      posteriorMethod = "numerical",
      prior = prior
    ),
    nterms = 4,
    nPiGrid = 31,
    nShapeGrid = 31
  )

  expect_identical(fit$method, "bayes")
  expect_identical(fit$posteriorMethod, "numerical")
  expect_s3_class(
    fit$posterior$representation,
    "numericalPosteriorRepresentation"
  )
  expect_s3_class(fit$posterior, "psPosterior")
})
