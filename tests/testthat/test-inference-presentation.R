test_that("predict preserves plug-in behaviour by default", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(20, 5, 1),
    type = "P"
  )
  fit = fitDist(pData, nterms = 4)

  expect_equal(predict(fit), predict(fit, type = "plugIn"))
  expect_equal(
    predict(fit, newdata = 0:2),
    fitted(fit, n = c(0, 1, 2), type = "plugIn")
  )
})

test_that("predict returns posterior means and credible intervals", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(8, 3, 1),
    type = "P"
  )
  fit = fitZIDist(
    pData,
    nterms = 3,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "numerical"),
    nPiGrid = 21,
    nShapeGrid = 21
  )

  expected = posteriorProbs(fit)
  predicted = predict(fit, type = "posteriorMean")
  interval = predict(
    fit,
    type = "posteriorMean",
    interval = "credible"
  )

  expect_equal(unname(predicted), expected$estimate)
  expect_identical(names(predicted), expected$term)
  expect_equal(interval$predicted, expected$estimate)
  expect_equal(interval$lower, expected$lower)
  expect_equal(interval$upper, expected$upper)
  expect_identical(rownames(interval), expected$term)
})

test_that("predict returns bootstrap means and percentile intervals", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(12, 6, 2),
    type = "P"
  )
  fit = fitZIDist(pData, nterms = 3)
  fit = bootstrapFit(
    fit,
    B = 30,
    seed = 931,
    silent = TRUE,
    parallel = FALSE
  )

  expected = bootstrapProbs(fit)
  predicted = predict(fit, type = "bootstrapMean")
  interval = predict(
    fit,
    type = "bootstrapMean",
    interval = "percentile"
  )

  expect_equal(unname(predicted), expected$estimate)
  expect_identical(names(predicted), expected$term)
  expect_equal(interval$predicted, expected$estimate)
  expect_equal(interval$lower, expected$lower)
  expect_equal(interval$upper, expected$upper)
  expect_identical(rownames(interval), expected$term)
})

test_that("predict selects stored uncertainty summaries by requested index", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(8, 3, 1),
    type = "P"
  )
  fit = fitZIDist(
    pData,
    nterms = 4,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "numerical"),
    nPiGrid = 21,
    nShapeGrid = 21
  )

  predicted = predict(
    fit,
    newdata = c(0, 2),
    type = "posteriorMean"
  )

  expect_identical(names(predicted), c("P0", "P2"))
  expect_error(
    predict(fit, newdata = 5, type = "posteriorMean"),
    "not available"
  )
})

test_that("predict rejects incompatible uncertainty intervals", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(8, 3, 1),
    type = "P"
  )
  fit = fitZIDist(
    pData,
    nterms = 3,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "numerical"),
    nPiGrid = 21,
    nShapeGrid = 21
  )

  expect_error(
    predict(fit, type = "posteriorMean", interval = "wald"),
    "posteriorMean"
  )
  expect_error(
    predict(fit, type = "plugIn", interval = "credible"),
    "plugIn"
  )
  expect_error(
    predict(fit, type = "posteriorMean", interval = "credible", level = 0.9),
    "stored at level"
  )
})

test_that("print psFit routes all fitted terms through fitted", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(12, 6, 2),
    type = "P"
  )
  fit = fitZIDist(pData, nterms = 3)

  expect_output(print(fit, nterms = 12), "P11")
})

test_that("psFit presentation reports attached bootstrap uncertainty", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(12, 6, 2),
    type = "P"
  )
  fit = fitZIDist(pData, nterms = 3)
  fit = bootstrapFit(
    fit,
    B = 25,
    seed = 932,
    silent = TRUE,
    parallel = FALSE
  )

  expect_output(print(fit), "Bootstrap uncertainty is attached")
  expect_output(summary(fit, nterms = 2), "Bootstrap probability summaries")
})

test_that("posterior and bootstrap summaries can be limited for presentation", {
  probabilitySummary = data.frame(
    term = paste0("P", 0:11),
    estimate = rep(0.05, 12),
    sd = rep(0.01, 12),
    lower = rep(0.03, 12),
    upper = rep(0.07, 12),
    level = rep(0.95, 12),
    posteriorMethod = rep("numerical", 12),
    stringsAsFactors = FALSE
  )
  posterior = newPsPosterior(
    method = "numerical",
    parameters = data.frame(
      parameter = c("pi", "shape"),
      estimate = c(0.2, 2),
      sd = c(0.05, 0.1)
    ),
    probabilities = probabilitySummary,
    representation = list(),
    level = 0.95
  )

  expect_output(print(posterior), "Showing 10 of 12")
  expect_equal(nrow(summary(posterior, nterms = 3)$probabilities), 3)

  bootstrapSummary = probabilitySummary
  bootstrapSummary$posteriorMethod = NULL
  bootstrapSummary$bootstrapMethod = "nonparametric"
  bootstrap = newPsBootstrap(
    method = "nonparametric",
    parameters = data.frame(
      parameter = c("pi", "shape"),
      estimate = c(0.2, 2),
      sd = c(0.05, 0.1),
      lower = c(0.1, 1.8),
      upper = c(0.3, 2.2),
      level = c(0.95, 0.95)
    ),
    probabilities = bootstrapSummary,
    replicates = data.frame(pi = 0.2, shape = 2),
    level = 0.95
  )

  expect_output(print(bootstrap), "Showing 10 of 12")
  expect_equal(nrow(summary(bootstrap, nterms = 3)$probabilities), 3)
})
