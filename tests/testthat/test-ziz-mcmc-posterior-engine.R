test_that("zero-inflated zeta MCMC is reproducible with a seed", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 4))

  fit1 = fitZIDist(
    pData,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "mcmc", prior = prior),
    shape1 = 2,
    shape2 = 5,
    theta0 = c(0.3, 2),
    nIter = 1200,
    nBurnIn = 300,
    seed = 779
  )
  fit2 = fitZIDist(
    pData,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "mcmc", prior = prior),
    shape1 = 2,
    shape2 = 5,
    theta0 = c(0.3, 2),
    nIter = 1200,
    nBurnIn = 300,
    seed = 779
  )

  expect_identical(fit1$chain, fit2$chain)
  expect_identical(fit1$fit$acceptance, fit2$fit$acceptance)
  expect_equal(fit1$pi, fit2$pi)
  expect_equal(fit1$shape, fit2$shape)
})

test_that("MCMC agrees with the numerical posterior under a non-uniform beta prior", {
  pData = makePSData(n = c(0, 1, 2, 3), count = c(12, 5, 2, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 4))

  numericalFit = fitZIDist(
    pData,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "numerical", prior = prior),
    shape1 = 2,
    shape2 = 5,
    nPiGrid = 81,
    nShapeGrid = 81
  )
  mcmcFit = fitZIDist(
    pData,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "mcmc", prior = prior),
    shape1 = 2,
    shape2 = 5,
    theta0 = c(0.3, 2),
    nIter = 12000,
    nBurnIn = 2000,
    seed = 20260806
  )

  expect_equal(mcmcFit$pi, numericalFit$pi, tolerance = 0.04)
  expect_equal(mcmcFit$shape, numericalFit$shape, tolerance = 0.12)
  piVarianceDifference = abs(
    unname(mcmcFit$var.cov["pi", "pi"]) -
      unname(numericalFit$var.cov["pi", "pi"])
  )
  expect_lte(piVarianceDifference, 0.001)
  expect_equal(
    unname(mcmcFit$var.cov["shape", "shape"]),
    unname(numericalFit$var.cov["shape", "shape"]),
    tolerance = 0.08
  )
  expect_true(all(mcmcFit$fit$acceptance > 0))
  expect_true(all(mcmcFit$fit$acceptance < 1))
})

test_that("legacy MCMC shape bounds define the prior actually used", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")

  fit = fitZIDist(
    pData,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "mcmc"),
    a = -1,
    b = 1,
    theta0 = c(0.4, 2),
    nIter = 1000,
    nBurnIn = 200,
    seed = 31
  )

  expect_equal(fit$bayesOptions$prior$range, 1 + exp(c(-1, 1)))
  expect_true(all(fit$chain$shape > fit$bayesOptions$prior$range[1]))
  expect_true(all(fit$chain$shape < fit$bayesOptions$prior$range[2]))
})

legacyMcmcZizFit = function(x,
                             prior,
                             nterms,
                             theta0 = c(0.5, 2),
                             shape1 = 1,
                             shape2 = 1,
                             nIter = 1000,
                             nBurnIn = 200,
                             seed = NULL,
                             level = 0.95) {
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  obsData = zizObservationData(x)
  counts = x$data$rn
  nTotal = nIter + nBurnIn
  currentPi = unname(theta0[1])
  currentShape = unname(theta0[2])
  chain = matrix(
    NA_real_,
    nrow = nIter,
    ncol = 2L,
    dimnames = list(NULL, c("pi", "shape"))
  )
  updateShape = sample(c(TRUE, FALSE), nTotal, replace = TRUE)
  logUniforms = log(runif(nTotal))
  acceptedShape = 0L
  acceptedPi = 0L
  proposedShape = 0L
  proposedPi = 0L

  for (iteration in seq_len(nTotal)) {
    if (updateShape[iteration]) {
      proposedShape = proposedShape + 1L
      candidateShape = runif(1L, prior$range[1], prior$range[2])

      logAcceptance =
        zizLogLikelihood(obsData, counts, currentPi, candidateShape) -
        zizLogLikelihood(obsData, counts, currentPi, currentShape) +
        prior$logd(candidateShape) - prior$logd(currentShape)

      if (is.finite(logAcceptance) &&
          (logAcceptance >= 0 || logUniforms[iteration] < logAcceptance)) {
        currentShape = candidateShape
        acceptedShape = acceptedShape + 1L
      }
    } else {
      proposedPi = proposedPi + 1L
      candidatePi = rbeta(1L, shape1 = shape1, shape2 = shape2)

      logTargetRatio =
        zizLogLikelihood(obsData, counts, candidatePi, currentShape) -
        zizLogLikelihood(obsData, counts, currentPi, currentShape) +
        dbeta(candidatePi, shape1, shape2, log = TRUE) -
        dbeta(currentPi, shape1, shape2, log = TRUE)
      logProposalRatio =
        dbeta(currentPi, shape1, shape2, log = TRUE) -
        dbeta(candidatePi, shape1, shape2, log = TRUE)
      logAcceptance = logTargetRatio + logProposalRatio

      if (is.finite(logAcceptance) &&
          (logAcceptance >= 0 || logUniforms[iteration] < logAcceptance)) {
        currentPi = candidatePi
        acceptedPi = acceptedPi + 1L
      }
    }

    if (iteration > nBurnIn) {
      chain[iteration - nBurnIn, ] = c(currentPi, currentShape)
    }
  }

  chain = as.data.frame(chain)
  par = c(pi = mean(chain$pi), shape = mean(chain$shape))
  posteriorProbs = summariseZizSampleProbabilities(
    pi = chain$pi,
    shape = chain$shape,
    type = x$type,
    nterms = nterms,
    level = level,
    posteriorMethod = "mcmc"
  )
  nvals = seq_len(nterms)
  fitted = (1 - par[["pi"]]) * dzetaStandard(nvals, shape = par[["shape"]])
  fitted[nvals == 1L] = fitted[nvals == 1L] + par[["pi"]]
  names(fitted) = if (x$type == "P") {
    paste0("P", nvals - 1L)
  } else {
    paste0("S", nvals)
  }

  list(
    par = par,
    acceptance = c(
      pi = if (proposedPi > 0L) acceptedPi / proposedPi else NA_real_,
      shape = if (proposedShape > 0L) acceptedShape / proposedShape else NA_real_
    ),
    var.cov = cov(chain),
    fitted = fitted,
    chain = chain,
    posteriorProbs = posteriorProbs
  )
}

test_that("MCMC ZIZ engine implements the posterior protocol", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 4))
  engine = mcmcPosteriorEngine()
  model = zizModel()

  representation = fitPosterior(
    engine = engine,
    model = model,
    x = pData,
    prior = prior,
    theta0 = c(0.3, 2),
    shape1 = 2,
    shape2 = 5,
    nIter = 1000,
    nBurnIn = 200,
    seed = 779,
    silent = TRUE
  )
  summary = summarisePosterior(engine, model, representation)
  pointEstimate = posteriorPointEstimate(engine, model, representation)
  diagnostics = posteriorDiagnostics(engine, representation)

  expect_s3_class(representation, "mcmcPosteriorRepresentation")
  expect_s3_class(representation$value$chain, "data.frame")
  expect_identical(names(representation$value$chain), c("pi", "shape"))
  expect_identical(summary$parameter, c("pi", "shape"))
  expect_equal(unname(pointEstimate), summary$estimate)
  expect_named(pointEstimate, c("pi", "shape"))
  expect_equal(representation$value$mean, pointEstimate, tolerance = 0)
  expect_equal(representation$value$variance, cov(representation$value$chain), tolerance = 0)
  expect_identical(diagnostics$model, "ziz")
  expect_identical(diagnostics$nIter, 1000L)
  expect_identical(diagnostics$nBurnIn, 200L)
  expect_true(all(diagnostics$acceptance > 0))
  expect_true(all(diagnostics$acceptance < 1))
})

test_that("MCMC ZIZ fitting reproduces the legacy sampler and orchestration", {
  prior = makePrior(family = "uniform", range = c(1.1, 4))

  for (type in c("P", "S")) {
    data = if (type == "P") {
      makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = type)
    } else {
      makePSData(n = c(1, 2, 3), count = c(8, 3, 1), type = type)
    }

    expected = legacyMcmcZizFit(
      x = data,
      prior = prior,
      nterms = 4,
      theta0 = c(0.3, 2),
      shape1 = 2,
      shape2 = 5,
      nIter = 1000,
      nBurnIn = 200,
      seed = 9127
    )
    actual = fitZIDistBayes(
      x = data,
      prior = prior,
      nterms = 4,
      theta0 = c(0.3, 2),
      shape1 = 2,
      shape2 = 5,
      nIter = 1000,
      nBurnIn = 200,
      seed = 9127,
      silent = TRUE
    )

    expect_identical(actual$chain, expected$chain)
    expect_equal(actual$fit$par, expected$par, tolerance = 0)
    expect_equal(actual$fit$acceptance, expected$acceptance, tolerance = 0)
    expect_equal(actual$var.cov, expected$var.cov, tolerance = 0)
    expect_equal(actual$fitted, expected$fitted, tolerance = 1e-14)
    expect_null(dim(actual$fitted))
    expect_equal(actual$posteriorProbs, expected$posteriorProbs, tolerance = 0)
    expect_s3_class(actual$posteriorRepresentation, "mcmcPosteriorRepresentation")
    expect_identical(actual$posterior$representation, actual$posteriorRepresentation)
    expect_equal(actual$posterior$diagnostics, actual$fit$acceptance, tolerance = 0)
  }
})

test_that("fitZIDist MCMC Bayesian path uses the migrated engine", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")
  prior = makePrior(family = "uniform", range = c(1.1, 4))

  fit = fitZIDist(
    pData,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "mcmc", prior = prior),
    shape1 = 2,
    shape2 = 5,
    theta0 = c(0.3, 2),
    nIter = 1000,
    nBurnIn = 200,
    seed = 778,
    silent = TRUE
  )

  expect_identical(fit$posteriorMethod, "mcmc")
  expect_s3_class(fit$posteriorRepresentation, "mcmcPosteriorRepresentation")
  expect_s3_class(fit$posterior, "psPosterior")
  expect_identical(fit$posterior$representation, fit$posteriorRepresentation)
})
