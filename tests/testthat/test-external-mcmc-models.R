test_that("external Poisson fits through generic MCMC Bayesian orchestration", {
  data = makePSData(
    n = c(0, 1, 2, 3),
    count = c(30, 12, 5, 1),
    type = "P"
  )
  prior = list(shape = 2, rate = 1)

  fitObject = fit(
    data,
    model = externalPoissonModel(),
    nterms = 5,
    method = "bayes",
    prior = prior,
    bayesOptions = list(posteriorMethod = "mcmc"),
    nIter = 300L,
    nBurnIn = 100L,
    proposalScale = 0.2,
    seed = 779
  )

  expect_s3_class(fitObject, "psFit")
  expect_s3_class(fitObject$modelObject, "externalPoissonModel")
  expect_identical(fitObject$method, "bayes")
  expect_identical(fitObject$posteriorMethod, "mcmc")
  expect_true(is.finite(fitObject$lambda))
  expect_gt(fitObject$lambda, 0)
  expect_identical(fitObject$posterior$diagnostics$generic, TRUE)
  expect_equal(nrow(fitObject$posterior$representation$value$chain), 300L)
  expect_true(all(is.finite(fitted(fitObject))))
  expect_true(is.finite(DIC(fitObject)))
})


test_that("external generic MCMC requires an explicit model-specific prior", {
  data = makePSData(
    n = c(0, 1, 2),
    count = c(20, 8, 2),
    type = "P"
  )

  expect_error(
    fit(
      data,
      model = externalPoissonModel(),
      method = "bayes",
      nIter = 50L,
      nBurnIn = 20L
    ),
    "require an explicit model-specific prior"
  )
})


test_that("unsupported external posterior engines fail deliberately", {
  data = makePSData(
    n = c(0, 1, 2),
    count = c(20, 8, 2),
    type = "P"
  )

  expect_error(
    fit(
      data,
      model = externalPoissonModel(),
      method = "bayes",
      prior = list(shape = 2, rate = 1),
      bayesOptions = list(posteriorMethod = "laplace")
    ),
    "mcmc"
  )
})


test_that("external Poisson-normal uses the generic transformed MCMC path", {
  data = makePSData(
    n = 0:4,
    count = c(150, 120, 55, 20, 5),
    type = "P"
  )
  prior = list(muMean = 0, muSd = 2, sigmaScale = 1)

  fitObject = fit(
    data,
    model = externalPoissonNormalModel(),
    nterms = 5,
    method = "bayes",
    prior = prior,
    bayesOptions = list(posteriorMethod = "mcmc"),
    nIter = 40L,
    nBurnIn = 15L,
    proposalScale = c(mu = 0.08, sigma = 0.08),
    seed = 779
  )

  expect_s3_class(fitObject$modelObject, "externalPoissonNormalModel")
  expect_identical(fitObject$posteriorMethod, "mcmc")
  expect_true(is.finite(fitObject$mu))
  expect_true(is.finite(fitObject$sigma))
  expect_gt(fitObject$sigma, 0)

  chain = fitObject$posterior$representation$value$chain
  expect_identical(names(chain), c("mu", "sigma"))
  expect_true(all(chain$sigma > 0))
})


test_that("external Poisson-normal Bayesian fits complete the public posterior contract", {
  data = makePSData(
    n = 0:4,
    count = c(150, 120, 55, 20, 5),
    type = "P"
  )
  prior = list(muMean = 0, muSd = 2, sigmaScale = 1)

  fitObject = fit(
    data,
    model = externalPoissonNormalModel(),
    nterms = 5,
    method = "bayes",
    prior = prior,
    bayesOptions = list(posteriorMethod = "mcmc"),
    nIter = 30L,
    nBurnIn = 10L,
    proposalScale = c(mu = 0.08, sigma = 0.08),
    seed = 779
  )

  posteriorSummary = summary(fitObject)
  expect_s3_class(posteriorSummary, "summary.psPosterior")
  expect_identical(
    posteriorSummary$parameters$parameter,
    c("mu", "sigma")
  )
  expect_true(all(is.finite(posteriorSummary$parameters$estimate)))
  expect_true(all(is.finite(posteriorSummary$parameters$sd)))

  probabilitySummary = fitObject$posterior$probabilities
  expect_identical(probabilitySummary$term, paste0("P", 0:4))
  expect_true(all(is.finite(probabilitySummary$estimate)))
  expect_true(all(probabilitySummary$estimate >= 0))
  expect_true(all(probabilitySummary$estimate <= 1))
  expect_true(all(probabilitySummary$lower <= probabilitySummary$estimate))
  expect_true(all(probabilitySummary$estimate <= probabilitySummary$upper))

  expectedFitted = externalPoissonNormalProbability(
    0:4,
    mu = fitObject$mu,
    sigma = fitObject$sigma
  )
  expect_equal(unname(fitted(fitObject)), expectedFitted, tolerance = 1e-8)
  expect_true(is.finite(DIC(fitObject)))

  path = tempfile(fileext = ".rds")
  saveRDS(fitObject, path)
  restored = readRDS(path)
  expect_s3_class(restored$modelObject, "externalPoissonNormalModel")
  expect_identical(restored$posteriorMethod, "mcmc")
  expect_equal(restored$posterior$parameters, fitObject$posterior$parameters)
  expect_equal(restored$posterior$probabilities, fitObject$posterior$probabilities)
  expect_equal(restored$posterior$representation$value$chain,
               fitObject$posterior$representation$value$chain)
  expect_equal(fitted(restored), fitted(fitObject))
})


test_that("external Poisson-normal Bayesian fitting preserves P and S support mapping", {
  counts = c(150, 120, 55, 20, 5)
  pData = makePSData(n = 0:4, count = counts, type = "P")
  sData = makePSData(n = 1:5, count = counts, type = "S")
  prior = list(muMean = 0, muSd = 2, sigmaScale = 1)

  fitExternal = function(data) {
    fit(
      data,
      model = externalPoissonNormalModel(),
      nterms = 5,
      method = "bayes",
      prior = prior,
      bayesOptions = list(posteriorMethod = "mcmc"),
      nIter = 25L,
      nBurnIn = 10L,
      proposalScale = c(mu = 0.08, sigma = 0.08),
      seed = 779
    )
  }

  pFit = fitExternal(pData)
  sFit = fitExternal(sData)

  expect_equal(pFit$mu, sFit$mu, tolerance = 1e-12)
  expect_equal(pFit$sigma, sFit$sigma, tolerance = 1e-12)
  expect_equal(
    pFit$posterior$representation$value$chain,
    sFit$posterior$representation$value$chain,
    tolerance = 1e-12
  )
  expect_equal(
    pFit$posterior$probabilities$estimate,
    sFit$posterior$probabilities$estimate,
    tolerance = 1e-12
  )
  expect_equal(unname(fitted(pFit)), unname(fitted(sFit)), tolerance = 1e-12)
})
