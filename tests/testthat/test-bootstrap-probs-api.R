test_that("fit method bootstrap attaches a public psBootstrap object", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(12, 6, 2),
    type = "P"
  )
  fitObject = fit(
    pData,
    model = zizModel(),
    method = "bootstrap",
    nterms = 3,
    B = 30,
    seed = 921,
    silent = TRUE,
    parallel = FALSE
  )

  expect_s3_class(fitObject, "psFit")
  expect_identical(fitObject$method, "bootstrap")
  expect_s3_class(fitObject$bootstrap, "psBootstrap")
  expect_identical(fitObject$bootstrap$diagnostics$B, 30L)
})

test_that("bootstrapProbs dispatches on psFit and psBootstrap", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(12, 6, 2),
    type = "P"
  )
  fitObject = fit(
    pData,
    model = zizModel(),
    method = "bootstrap",
    nterms = 3,
    B = 35,
    seed = 922,
    silent = TRUE,
    parallel = FALSE
  )

  expect_equal(
    bootstrapProbs(fitObject),
    bootstrapProbs(fitObject$bootstrap)
  )
  expect_equal(bootstrapProbs(fitObject, n = 2)$term, c("P0", "P1"))
  expect_equal(bootstrapProbs(fitObject, n = c(0, 2))$term, c("P0", "P2"))
})

test_that("bootstrapProbs uses S indexing rules", {
  sData = makePSData(
    n = c(1, 2, 3),
    count = c(12, 5, 2),
    type = "S"
  )
  fitObject = fit(
    sData,
    model = zizModel(),
    method = "bootstrap",
    nterms = 3,
    B = 30,
    seed = 923,
    silent = TRUE,
    parallel = FALSE
  )

  expect_equal(bootstrapProbs(fitObject, n = c(1, 3))$term, c("S1", "S3"))
  expect_error(bootstrapProbs(fitObject, n = c(0, 1)), "S indices")
})

test_that("bootstrapProbs requires bootstrap inference", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(12, 6, 2),
    type = "P"
  )
  fitObject = fit(pData, model = zizModel(), nterms = 3)

  expect_error(bootstrapProbs(fitObject), "method = 'bootstrap'")
  expect_error(
    fitted(fitObject, type = "bootstrapMean"),
    "method = 'bootstrap'"
  )
})

test_that("fitted bootstrapMean returns bootstrap probability means", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(12, 6, 2),
    type = "P"
  )
  mleFit = fit(pData, model = zizModel(), nterms = 3)
  plugIn = fitted(mleFit)
  fitObject = fit(
    pData,
    model = zizModel(),
    method = "bootstrap",
    nterms = 3,
    B = 40,
    seed = 924,
    silent = TRUE,
    parallel = FALSE
  )

  expected = bootstrapProbs(fitObject)$estimate
  names(expected) = bootstrapProbs(fitObject)$term

  expect_equal(fitted(fitObject, type = "bootstrapMean"), expected)
  expect_equal(fitted(fitObject), plugIn)
  expect_equal(fitted(fitObject, type = "plugIn"), plugIn)
})

test_that("plot.psBootstrap returns plotted summaries invisibly", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(12, 6, 2),
    type = "P"
  )
  fitObject = fit(
    pData,
    model = zizModel(),
    method = "bootstrap",
    nterms = 3,
    B = 25,
    seed = 925,
    silent = TRUE,
    parallel = FALSE
  )

  plotFile = tempfile(fileext = ".pdf")
  grDevices::pdf(plotFile)
  plotted = plot(fitObject$bootstrap, n = 2)
  grDevices::dev.off()

  expect_equal(plotted$term, c("P0", "P1"))
  expect_true(file.exists(plotFile))
})

test_that("fit method bootstrap reports unsupported model coverage", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(12, 6, 2),
    type = "P"
  )

  expect_error(
    fit(pData, model = logarithmicModel(), method = "bootstrap", B = 10),
    "currently supports zeta and ziz"
  )
})

test_that("bootstrapFit retains its historical MLE-fit requirement", {
  oldOption = getOption("fitPS.deprecationWarnings")
  on.exit(options(fitPS.deprecationWarnings = oldOption), add = TRUE)
  options(fitPS.deprecationWarnings = FALSE)

  pData = makePSData(
    n = c(0, 1, 2),
    count = c(12, 6, 2),
    type = "P"
  )
  bayesFit = structure(
    list(
      psData = pData,
      method = "bayes",
      model = "ziz",
      fitted = c(P0 = 0.4, P1 = 0.3, P2 = 0.2)
    ),
    class = "psFit"
  )

  expect_error(bootstrapFit(bayesFit, B = 10), "MLE psFit")
})
