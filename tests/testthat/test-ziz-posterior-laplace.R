test_that("zero-inflated Laplace approximation returns finite posterior mode", {
  pData = makePSData(n = c(0, 1, 2), count = c(7, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 5))

  approximation = makeZizPosteriorLaplace(
    x = pData,
    prior = prior,
    start = c(pi = 0.4, shape = 2)
  )

  expect_named(approximation$mode, c("pi", "shape"))
  expect_named(approximation$modeWorking, c("eta", "tau"))
  expect_true(approximation$mode[["pi"]] > 0)
  expect_true(approximation$mode[["pi"]] < 1)
  expect_true(approximation$mode[["shape"]] > 1)
  expect_true(is.finite(approximation$logPosteriorMode))
  expect_equal(dim(approximation$covarianceWorking), c(2L, 2L))
  expect_equal(dim(approximation$varCov), c(2L, 2L))
})

test_that("zero-inflated Laplace mode improves on the numerical grid mean", {
  pData = makePSData(n = c(0, 1, 2), count = c(7, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 5))

  grid = makeZizPosteriorGrid(
    x = pData,
    prior = prior,
    nPiGrid = 31,
    nShapeGrid = 31
  )
  approximation = makeZizPosteriorLaplace(
    x = pData,
    prior = prior,
    start = grid$mean
  )

  obsData = zizObservationData(pData)
  counts = pData$data$rn
  gridMeanWorking = zizThetaToWorking(grid$mean)
  gridMeanLogPosterior = zizWorkingLogPosterior(
    working = gridMeanWorking,
    obsData = obsData,
    counts = counts,
    prior = prior
  )

  expect_true(is.finite(gridMeanLogPosterior))
  expect_true(approximation$logPosteriorMode >= gridMeanLogPosterior - 1e-8)
  expect_true(approximation$mode[["pi"]] > min(grid$pi))
  expect_true(approximation$mode[["pi"]] < max(grid$pi))
  expect_true(approximation$mode[["shape"]] >= min(grid$shape))
  expect_true(approximation$mode[["shape"]] <= max(grid$shape))
})

test_that("Laplace posterior engine preserves ZIZ fit orchestration", {
  prior = makePrior(family = "uniform", range = c(1.1, 5))

  for (type in c("P", "S")) {
    x = if (type == "P") {
      makePSData(n = c(0, 1, 2), count = c(7, 3, 1), type = "P")
    } else {
      makePSData(n = c(1, 2, 3), count = c(7, 3, 1), type = "S")
    }

    approximation = makeZizPosteriorLaplace(
      x = x,
      prior = prior,
      start = c(pi = 0.4, shape = 2)
    )
    set.seed(314)
    workingDraws = makeZizProposalDraws(
      mean = approximation$modeWorking,
      covariance = approximation$covarianceWorking,
      n = 300
    )
    thetaDraws = t(apply(workingDraws, 1L, zizWorkingToTheta))
    expectedPosteriorProbs = summariseZizSampleProbabilities(
      pi = thetaDraws[, "pi"],
      shape = thetaDraws[, "shape"],
      type = x$type,
      nterms = 4,
      level = 0.95,
      posteriorMethod = "laplace"
    )
    nvals = 1:4
    expectedFitted = (1 - approximation$mode[["pi"]]) *
      dzetaStandard(nvals, shape = approximation$mode[["shape"]])
    expectedFitted[nvals == 1] = expectedFitted[nvals == 1] +
      approximation$mode[["pi"]]
    names(expectedFitted) = if (type == "P") {
      paste0("P", nvals - 1)
    } else {
      paste0("S", nvals)
    }

    actual = fitZIDistBayesLaplace(
      x = x,
      nterms = 4,
      prior = prior,
      start = c(pi = 0.4, shape = 2),
      nPosteriorDraws = 300,
      seed = 314
    )

    expect_s3_class(actual$posteriorRepresentation, "laplacePosteriorRepresentation")
    expect_equal(actual$fit$par, approximation$mode)
    expect_equal(actual$var.cov, approximation$varCov)
    expect_equal(actual$fitted, expectedFitted)
    expect_null(dim(actual$fitted))
    expect_equal(actual$posteriorProbs, expectedPosteriorProbs)
    expect_equal(actual$laplace, approximation)
  }
})

test_that("Laplace engine is available only for supported models", {
  pData = makePSData(n = c(0, 1, 2), count = c(7, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 5))

  representation = fitPosterior(
    laplacePosteriorEngine(),
    zizModel(),
    pData,
    prior,
    start = c(pi = 0.4, shape = 2)
  )

  expect_s3_class(representation, "laplacePosteriorRepresentation")
  expect_named(
    posteriorPointEstimate(laplacePosteriorEngine(), zizModel(), representation),
    c("pi", "shape")
  )
  expect_error(
    fitPosterior(laplacePosteriorEngine(), zetaModel(), pData, prior),
    "not supported"
  )
})
