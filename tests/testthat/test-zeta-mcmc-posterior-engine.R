legacyMcmcZetaFit = function(x,
                              prior,
                              nterms,
                              shape0 = 2,
                              nIter = 1000,
                              nBurnIn = 100,
                              silent = TRUE) {
  nvals = seq_len(nterms)
  observations = if (x$type == "P") {
    x$data$n + 1
  } else {
    x$data$n
  }

  logLikelihood = function(shape) {
    sum(x$data$rn * dzetaStandard(observations, shape = shape, log = TRUE))
  }

  nTotal = nIter + nBurnIn
  draws = runif(nTotal, prior$range[1], prior$range[2])
  chain = numeric(nIter)
  logUniforms = log(runif(nTotal))

  currentLogPosterior = logLikelihood(shape0) + prior$logd(shape0)
  if (!is.finite(currentLogPosterior)) {
    stop("Log likelihood is not finite at starting value")
  }

  i = 1
  while (i <= nTotal) {
    proposedShape = draws[i]
    proposedLogPosterior = logLikelihood(proposedShape) + prior$logd(proposedShape)

    if (proposedLogPosterior > currentLogPosterior ||
        logUniforms[i] < (proposedLogPosterior - currentLogPosterior)) {
      shape0 = proposedShape
      currentLogPosterior = proposedLogPosterior
    }

    if (i > nBurnIn) {
      chain[i - nBurnIn] = shape0
    }
    i = i + 1
  }

  shape = mean(chain)
  variance = var(chain)
  fitted = dzetaStandard(nvals, shape = shape)
  names(fitted) = if (x$type == "P") {
    paste0("P", nvals - 1)
  } else {
    paste0("S", nvals)
  }

  densityEstimate = density(chain, from = prior$range[1])

  list(
    shape = shape,
    var.shape = variance,
    fitted = fitted,
    chain = chain,
    pdf = splinefun(densityEstimate$x, densityEstimate$y)
  )
}

test_that("MCMC zeta engine implements the posterior protocol", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 4))
  engine = mcmcPosteriorEngine()
  model = zetaModel()

  set.seed(42)
  representation = fitPosterior(
    engine,
    model,
    pData,
    prior,
    nIter = 1000,
    nBurnIn = 100,
    silent = TRUE
  )
  summary = summarisePosterior(engine, model, representation)
  pointEstimate = posteriorPointEstimate(engine, model, representation)
  diagnostics = posteriorDiagnostics(engine, representation)

  expect_s3_class(representation, "mcmcPosteriorRepresentation")
  expect_length(representation$value$chain, 1000)
  expect_named(summary, c("parameter", "estimate", "sd"))
  expect_identical(summary$parameter, "shape")
  expect_equal(unname(pointEstimate), summary$estimate)
  expect_named(pointEstimate, "shape")
  expect_true(summary$sd > 0)
  expect_identical(diagnostics$model, "zeta")
  expect_identical(diagnostics$nIter, 1000)
  expect_identical(diagnostics$nBurnIn, 100)
})

test_that("MCMC zeta fitting reproduces the legacy algorithm", {
  prior = makePrior(family = "uniform", range = c(1.1, 4))

  for (type in c("P", "S")) {
    data = if (type == "P") {
      makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = type)
    } else {
      makePSData(n = c(1, 2, 3), count = c(8, 3, 1), type = type)
    }

    set.seed(779)
    expected = legacyMcmcZetaFit(
      data,
      prior,
      nterms = 4,
      nIter = 1000,
      nBurnIn = 100
    )
    set.seed(779)
    actual = fitDist(
      data,
      prior = prior,
      nterms = 4,
      method = "bayes",
      bayesOptions = list(posteriorMethod = "mcmc"),
      nIter = 1000,
      nBurnIn = 100,
      silent = TRUE
    )

    expect_identical(actual$chain, expected$chain)
    expect_equal(actual$shape, expected$shape, tolerance = 0)
    expect_equal(actual$var.shape, expected$var.shape, tolerance = 1e-14)
    expect_equal(actual$fitted, expected$fitted, tolerance = 1e-14)
    grid = seq(1.2, 3.8, length.out = 7)
    expect_equal(actual$pdf(grid), expected$pdf(grid), tolerance = 1e-12)
    expect_s3_class(actual$posterior$representation, "mcmcPosteriorRepresentation")
  }
})

test_that("fitDist MCMC Bayesian path uses the migrated engine", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 4))

  set.seed(123)
  fit = fitDist(
    pData,
    method = "bayes",
    bayesOptions = list(
      posteriorMethod = "mcmc",
      prior = prior
    ),
    nterms = 4,
    nIter = 1000,
    nBurnIn = 100,
    silent = TRUE
  )

  expect_identical(fit$method, "bayes")
  expect_identical(fit$posteriorMethod, "mcmc")
  expect_s3_class(fit$posterior$representation, "mcmcPosteriorRepresentation")
  expect_length(fit$chain, 1000)
  expect_true(is.function(fit$pdf))
})

