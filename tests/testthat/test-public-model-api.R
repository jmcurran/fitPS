test_that("public model fits support common fitted-object presentation", {
  pData = makePSData(
    n = c(0, 1, 2, 3),
    count = c(30, 12, 5, 1),
    type = "P"
  )
  object = fit(pData, model = externalPoissonModel(), nterms = 5)

  summaryMatrix = NULL
  summaryOutput = capture.output({
    summaryMatrix = summary(object)
  })
  printOutput = capture.output(print(object, nterms = 3))

  expect_identical(rownames(summaryMatrix), "lambda")
  expect_identical(colnames(summaryMatrix), c("Estimate", "Std.Err"))
  expect_equal(summaryMatrix["lambda", "Estimate"], object$lambda)
  expect_true(is.finite(summaryMatrix["lambda", "Std.Err"]))
  expect_true(any(grepl("lambda", summaryOutput, fixed = TRUE)))
  expect_true(any(grepl("Estimated model parameters", printOutput, fixed = TRUE)))
  expect_true(any(grepl("fitted values", printOutput, fixed = TRUE)))
})


test_that("multi-parameter external fits use the same generic summary contract", {
  pData = makePSData(
    n = 0:8,
    count = c(3032, 3240, 2035, 997, 426, 168, 64, 24, 9),
    type = "P"
  )
  object = fit(pData, model = externalPoissonNormalModel(), nterms = 9)

  summaryMatrix = NULL
  capture.output({
    summaryMatrix = summary(object)
  })

  expect_identical(rownames(summaryMatrix), c("mu", "sigma"))
  expect_identical(colnames(summaryMatrix), c("Estimate", "Std.Err"))
  expect_equal(unname(summaryMatrix[, "Estimate"]), c(object$mu, object$sigma))
  expect_true(all(is.finite(summaryMatrix[, "Std.Err"]) | is.na(summaryMatrix[, "Std.Err"])))
})


test_that("external models reject Bayesian fitting until they declare a supported path", {
  pData = makePSData(n = 0:2, count = c(20, 5, 1), type = "P")

  expect_error(
    fit(pData, model = externalPoissonModel(), method = "bayes"),
    "generic external-model fitting currently supports method = 'mle'"
  )
})


test_that("external plug-in probability functions respect P and S support mapping", {
  counts = c(30, 12, 5, 1)
  pData = makePSData(n = 0:3, count = counts, type = "P")
  sData = makePSData(n = 1:4, count = counts, type = "S")

  pFit = fit(pData, model = externalPoissonModel(), nterms = 5)
  sFit = fit(sData, model = externalPoissonModel(), nterms = 5)

  expect_equal(
    unname(probfun(pFit)(0:4)),
    unname(probfun(sFit)(1:5)),
    tolerance = 1e-8
  )
})
