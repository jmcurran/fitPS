test_that("bootstrapFit attaches a public psBootstrap object", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(12, 6, 2),
    type = "P"
  )
  fit = fitZIDist(pData, nterms = 3)
  fit = bootstrapFit(
    fit,
    B = 30,
    seed = 921,
    silent = TRUE,
    parallel = FALSE
  )

  expect_s3_class(fit, "psFit")
  expect_s3_class(fit$bootstrap, "psBootstrap")
  expect_identical(fit$bootstrap$diagnostics$B, 30L)
})

test_that("bootstrapProbs dispatches on psFit and psBootstrap", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(12, 6, 2),
    type = "P"
  )
  fit = fitZIDist(pData, nterms = 3)
  fit = bootstrapFit(
    fit,
    B = 35,
    seed = 922,
    silent = TRUE,
    parallel = FALSE
  )

  expect_equal(bootstrapProbs(fit), bootstrapProbs(fit$bootstrap))
  expect_equal(bootstrapProbs(fit, n = 2)$term, c("P0", "P1"))
  expect_equal(bootstrapProbs(fit, n = c(0, 2))$term, c("P0", "P2"))
})

test_that("bootstrapProbs uses S indexing rules", {
  sData = makePSData(
    n = c(1, 2, 3),
    count = c(12, 5, 2),
    type = "S"
  )
  fit = fitZIDist(sData, nterms = 3)
  fit = bootstrapFit(
    fit,
    B = 30,
    seed = 923,
    silent = TRUE,
    parallel = FALSE
  )

  expect_equal(bootstrapProbs(fit, n = c(1, 3))$term, c("S1", "S3"))
  expect_error(bootstrapProbs(fit, n = c(0, 1)), "S indices")
})

test_that("bootstrapProbs requires an attached bootstrap", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(12, 6, 2),
    type = "P"
  )
  fit = fitZIDist(pData, nterms = 3)

  expect_error(bootstrapProbs(fit), "bootstrapFit")
  expect_error(fitted(fit, type = "bootstrapMean"), "bootstrapFit")
})

test_that("fitted bootstrapMean returns bootstrap probability means", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(12, 6, 2),
    type = "P"
  )
  fit = fitZIDist(pData, nterms = 3)
  plugIn = fitted(fit)
  fit = bootstrapFit(
    fit,
    B = 40,
    seed = 924,
    silent = TRUE,
    parallel = FALSE
  )

  expected = bootstrapProbs(fit)$estimate
  names(expected) = bootstrapProbs(fit)$term

  expect_equal(fitted(fit, type = "bootstrapMean"), expected)
  expect_equal(fitted(fit), plugIn)
  expect_equal(fitted(fit, type = "plugIn"), plugIn)
})

test_that("plot.psBootstrap returns plotted summaries invisibly", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(12, 6, 2),
    type = "P"
  )
  fit = fitZIDist(pData, nterms = 3)
  fit = bootstrapFit(
    fit,
    B = 25,
    seed = 925,
    silent = TRUE,
    parallel = FALSE
  )

  plotFile = tempfile(fileext = ".pdf")
  grDevices::pdf(plotFile)
  plotted = plot(fit$bootstrap, n = 2)
  grDevices::dev.off()

  expect_equal(plotted$term, c("P0", "P1"))
  expect_true(file.exists(plotFile))
})

test_that("bootstrapFit rejects Bayesian fits", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(12, 6, 2),
    type = "P"
  )
  fit = structure(
    list(
      psData = pData,
      method = "bayes",
      model = "ziz",
      fitted = c(P0 = 0.4, P1 = 0.3, P2 = 0.2)
    ),
    class = "psFit"
  )

  expect_error(bootstrapFit(fit, B = 10), "MLE psFit")
})
