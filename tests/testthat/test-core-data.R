test_that("makePSData creates P data from counts", {
  psData = makePSData(n = 0:2, count = c(8, 3, 1), type = "P", notes = "source")

  expect_s3_class(psData, "psData")
  expect_equal(psData$type, "P")
  expect_equal(psData$data, data.frame(n = 0:2, rn = c(8, 3, 1)))
  expect_equal(psData$notes, "source")
})

test_that("makePSData tabulates raw observations", {
  psData = makePSData(n = c(0, 0, 1, 2, 2, 2), type = "P")

  expect_s3_class(psData, "psData")
  expect_equal(psData$data, data.frame(n = 0:2, rn = c(2, 1, 3)))
})

test_that("readData reads a temporary P csv", {
  csvFile = tempfile(fileext = ".csv")
  write.csv(data.frame(P = 0:2, count = c(8, 3, 1)), csvFile, row.names = FALSE)

  psData = readData(csvFile, notes = "temporary")

  expect_s3_class(psData, "psData")
  expect_equal(psData$type, "P")
  expect_equal(psData$data, data.frame(n = 0:2, rn = c(8L, 3L, 1L)))
  expect_equal(psData$notes, "temporary")
})

test_that("fitDist fits tiny simple P data", {
  psData = makePSData(n = 0:2, count = c(8, 3, 1), type = "P")

  fit = fitDist(psData, nterms = 3)

  expect_s3_class(fit, "psFit")
  expect_equal(fit$model, "zeta")
  expect_equal(fit$method, "mle")
  expect_true(is.finite(fit$shape))
  expect_named(fit$fitted, c("P0", "P1", "P2"))
  expect_true(all(fit$fitted >= 0))
})

test_that("predict.psFit returns fitted and new zeta predictions", {
  psData = makePSData(n = 0:2, count = c(8, 3, 1), type = "P")
  fit = fitDist(psData, nterms = 3)

  expect_equal(predict(fit), fit$fitted)

  predictions = predict(fit, newdata = 0:2)

  expect_named(predictions, c("P0", "P1", "P2"))
  expect_equal(unname(predictions), unname(fit$fitted))
})
