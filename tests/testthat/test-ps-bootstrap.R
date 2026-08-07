test_that("bootstrapPsFit creates a reproducible zeta psBootstrap object", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(10, 5, 2),
    type = "P"
  )
  fit = fitDist(pData, nterms = 3)

  first = bootstrapPsFit(
    fit,
    B = 40,
    seed = 812,
    silent = TRUE,
    parallel = FALSE
  )
  second = bootstrapPsFit(
    fit,
    B = 40,
    seed = 812,
    silent = TRUE,
    parallel = FALSE
  )

  expect_s3_class(first$bootstrap, "psBootstrap")
  expect_equal(first$bootstrap$replicates, second$bootstrap$replicates)
  expect_equal(first$bootstrap$probabilities, second$bootstrap$probabilities)
  expect_equal(first$bootstrap$probabilities$term, c("P0", "P1", "P2"))
  expect_identical(first$bootstrap$method, "nonparametric")
})

test_that("zero-inflated bootstrap summaries transform successful replicates", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(12, 6, 2),
    type = "P"
  )
  fit = fitZIDist(pData, nterms = 3)
  fit = bootstrapPsFit(
    fit,
    B = 50,
    seed = 813,
    silent = TRUE,
    parallel = FALSE
  )

  successful = complete.cases(fit$bootstrap$replicates)
  replicates = fit$bootstrap$replicates[successful, , drop = FALSE]
  transformed = zizProbabilities(
    pi = replicates$pi,
    shape = replicates$shape,
    n = 0:2,
    type = "P"
  )

  expect_equal(
    unname(fit$bootstrap$probabilities$estimate),
    unname(colMeans(transformed))
  )
  expect_equal(
    fit$bootstrap$diagnostics$nSuccessful +
      fit$bootstrap$diagnostics$nFailed,
    50
  )
  expect_true(all(fit$bootstrap$probabilities$lower >= 0))
  expect_true(all(fit$bootstrap$probabilities$upper <= 1))
})

test_that("bootstrap diagnostics retain failed zero-inflated resamples", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(1, 1, 1),
    type = "P"
  )
  fit = structure(
    list(
      psData = pData,
      fitted = c(P0 = 0.4, P1 = 0.35, P2 = 0.25),
      model = "ziz",
      method = "mle"
    ),
    class = "psFit"
  )
  fit = bootstrapPsFit(
    fit,
    B = 60,
    seed = 814,
    silent = TRUE,
    parallel = FALSE
  )

  expect_equal(nrow(fit$bootstrap$replicates), 60)
  expect_gt(fit$bootstrap$diagnostics$nFailed, 0)
  expect_equal(
    fit$bootstrap$diagnostics$failureRate,
    fit$bootstrap$diagnostics$nFailed / 60
  )
  expect_true(any(!complete.cases(fit$bootstrap$replicates)))
})

test_that("S bootstrap probabilities use valid S term names", {
  sData = makePSData(
    n = c(1, 2, 3),
    count = c(10, 4, 1),
    type = "S"
  )
  fit = fitZIDist(sData, nterms = 3)
  fit = bootstrapPsFit(
    fit,
    B = 40,
    seed = 815,
    silent = TRUE,
    parallel = FALSE
  )

  expect_equal(fit$bootstrap$probabilities$term, c("S1", "S2", "S3"))
  expect_equal(
    fit$bootstrap$parameters$parameter,
    c("pi", "shape")
  )
})

test_that("legacy bootFit retains its established return shapes", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(10, 5, 2),
    type = "P"
  )

  zetaValues = bootFit(
    pData,
    B = 20,
    model = "zeta",
    seed = 816,
    silent = TRUE,
    parallel = FALSE
  )
  zizValues = bootFit(
    pData,
    B = 20,
    model = "ziz",
    seed = 816,
    silent = TRUE,
    parallel = FALSE
  )

  expect_type(zetaValues, "double")
  expect_true(length(zetaValues) <= 20)
  expect_s3_class(zizValues, "data.frame")
  expect_named(zizValues, c("pi", "shape"))
  expect_true(nrow(zizValues) <= 20)
})

test_that("psBootstrap print and summary methods expose diagnostics", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(10, 5, 2),
    type = "P"
  )
  fit = fitDist(pData, nterms = 2)
  fit = bootstrapPsFit(
    fit,
    B = 20,
    seed = 817,
    silent = TRUE,
    parallel = FALSE
  )

  expect_output(print(fit$bootstrap), "fitPS bootstrap distribution")
  bootstrapSummary = summary(fit$bootstrap)
  expect_s3_class(bootstrapSummary, "summary.psBootstrap")
  expect_output(print(bootstrapSummary), "Summary of fitPS bootstrap")
})
