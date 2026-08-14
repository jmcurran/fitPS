test_that("weighted probability summaries retain weights", {
  probabilities = matrix(
    c(0.1, 0.9, 0.2, 0.8),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(NULL, c("P0", "P1"))
  )

  summary = summariseZizProbabilities(
    probabilities = probabilities,
    weights = c(0.75, 0.25),
    level = 0.8,
    posteriorMethod = "importance"
  )

  expect_equal(summary$estimate, c(0.125, 0.875))
  expect_equal(summary$posteriorMethod, rep("importance", 2L))
  expect_true(all(summary$lower >= 0 & summary$upper <= 1))
})

test_that("numerical posterior probabilities use grid weights", {
  data(Psurveys)

  fit = fitZIDist(
    Psurveys$roux,
    nterms = 4,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "numerical"),
    nPiGrid = 31,
    nShapeGrid = 31,
    shape1 = 2,
    shape2 = 5
  )

  expect_equal(fit$posterior$probabilities$term, paste0("P", 0:3))
  expect_equal(fit$posterior$probabilities$posteriorMethod, rep("numerical", 4L))
  expect_true(all(fit$posterior$probabilities$lower >= 0))
  expect_true(all(fit$posterior$probabilities$upper <= 1))
  expect_true(all(
    fit$posterior$probabilities$estimate >= fit$posterior$probabilities$lower &
      fit$posterior$probabilities$estimate <= fit$posterior$probabilities$upper
  ))
})

test_that("MCMC posterior probabilities transform retained draws", {
  data(Psurveys)

  fit = fitZIDist(
    Psurveys$roux,
    nterms = 3,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "mcmc"),
    nIter = 1000,
    nBurnIn = 300,
    seed = 20260806,
    shape1 = 2,
    shape2 = 5
  )

  transformed = zizProbabilities(
    pi = fit$posterior$representation$value$chain$pi,
    shape = fit$posterior$representation$value$chain$shape,
    n = 0:2,
    type = "P"
  )

  expect_equal(unname(fit$posterior$probabilities$estimate), unname(colMeans(transformed)))
  expect_equal(fit$posterior$probabilities$term, paste0("P", 0:2))
})

test_that("importance posterior probabilities use retained weights", {
  sData = makePSData(
    n = c(1, 2, 3),
    count = c(12, 4, 1),
    type = "S"
  )

  fit = fitZIDist(
    sData,
    nterms = 3,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "importance"),
    nSamples = 400,
    seed = 1234
  )

  transformed = zizProbabilities(
    pi = fit$posterior$representation$value$approximation$samples$pi,
    shape = fit$posterior$representation$value$approximation$samples$shape,
    n = 1:3,
    type = "S"
  )
  expected = colSums(transformed * fit$posterior$representation$value$approximation$samples$weight)

  expect_equal(unname(fit$posterior$probabilities$estimate), unname(expected))
  expect_equal(fit$posterior$probabilities$term, paste0("S", 1:3))
})

test_that("Laplace posterior probability summaries are seeded", {
  data(Psurveys)

  first = fitZIDist(
    Psurveys$roux,
    nterms = 3,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "laplace"),
    nPosteriorDraws = 500,
    seed = 99
  )
  second = fitZIDist(
    Psurveys$roux,
    nterms = 3,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "laplace"),
    nPosteriorDraws = 500,
    seed = 99
  )

  expect_identical(first$posterior, second$posterior)
  expect_equal(first$posterior$probabilities$posteriorMethod, rep("laplace", 3L))
})

test_that("posterior means can differ from plug-in probabilities", {
  pi = c(0.05, 0.95)
  shape = c(1.2, 4.5)
  transformed = zizProbabilities(pi, shape, n = 0:2, type = "P")
  posteriorMean = colMeans(transformed)
  plugIn = zizProbabilities(
    pi = mean(pi),
    shape = mean(shape),
    n = 0:2,
    type = "P"
  )[1L, ]

  expect_true(any(abs(posteriorMean - plugIn) > 0.01))
})
