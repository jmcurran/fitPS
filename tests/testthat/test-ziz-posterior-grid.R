test_that("zero-inflated zeta numerical posterior grid normalizes marginals", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 4))

  posteriorGrid = makeZizPosteriorGrid(
    pData,
    prior = prior,
    nPiGrid = 21,
    nShapeGrid = 21
  )

  expect_named(posteriorGrid$mean, c("pi", "shape"))
  expect_equal(
    sum(posteriorGrid$marginalDensity$pi) * posteriorGrid$dPi,
    1,
    tolerance = 0.08
  )
  expect_equal(
    sum(posteriorGrid$marginalDensity$shape) * posteriorGrid$dShape,
    1,
    tolerance = 0.08
  )
  expect_true(is.matrix(posteriorGrid$varCov))
  expect_named(as.data.frame(posteriorGrid$varCov), c("pi", "shape"))
})

test_that("fitZIDist uses numerical posterior as Bayesian default", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 4))

  fit = fitZIDist(
    pData,
    method = "bayes",
    bayesOptions = list(prior = prior),
    nPiGrid = 21,
    nShapeGrid = 21
  )

  expect_s3_class(fit, "psFit")
  expect_identical(fit$method, "bayes")
  expect_identical(fit$posteriorMethod, "numerical")
  expect_true(is.list(fit$posteriorGrid))
  expect_true(is.function(fit$marginalPdf$pi))
  expect_true(is.function(fit$marginalPdf$shape))
  expect_true(fit$pi > 0 && fit$pi < 1)
  expect_true(fit$shape > 1)
})

test_that("plotPosterior can use stored zero-inflated numerical marginals", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 4))
  fit = fitZIDist(
    pData,
    method = "bayes",
    bayesOptions = list(prior = prior),
    nPiGrid = 21,
    nShapeGrid = 21
  )

  pdfGrid = plotPosterior(fit, parameter = "pi", nGrid = 64)

  expect_s3_class(pdfGrid, "data.frame")
  expect_named(pdfGrid, c("x", "density"))
  expect_true(all(is.finite(pdfGrid$density)))
})
