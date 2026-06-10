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
