test_that("plot.psPosterior returns plotted probability summaries", {
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

  plotFile = tempfile(fileext = ".pdf")
  grDevices::pdf(plotFile)
  on.exit(grDevices::dev.off(), add = TRUE)

  plotted = plot(fit$posterior, n = 3)

  expect_s3_class(plotted, "data.frame")
  expect_equal(plotted$term, c("P0", "P1", "P2"))
  expect_equal(plotted, posteriorProbs(fit, n = 3))
})

test_that("plot.psPosterior can omit credible intervals", {
  posterior = newPsPosterior(
    method = "numerical",
    parameters = data.frame(
      parameter = c("pi", "shape"),
      estimate = c(0.2, 2.1),
      sd = c(0.05, 0.2)
    ),
    probabilities = data.frame(
      term = c("S1", "S2"),
      estimate = c(0.7, 0.2),
      sd = c(0.04, 0.03),
      lower = c(0.62, 0.14),
      upper = c(0.78, 0.26),
      level = c(0.95, 0.95),
      posteriorMethod = c("numerical", "numerical")
    ),
    representation = NULL
  )

  plotFile = tempfile(fileext = ".pdf")
  grDevices::pdf(plotFile)
  on.exit(grDevices::dev.off(), add = TRUE)

  plotted = plot(posterior, showInterval = FALSE)

  expect_equal(plotted$term, c("S1", "S2"))
})

test_that("plot.psPosterior validates its input", {
  expect_error(
    plot.psPosterior(list()),
    "class psPosterior"
  )
})
