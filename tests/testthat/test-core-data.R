test_that("makePSData constructs P and S survey data", {
  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")
  sData = makePSData(n = c(1, 2, 3), count = c(8, 3, 1), type = "S")

  expect_s3_class(pData, "psData")
  expect_s3_class(sData, "psData")
  expect_identical(pData$type, "P")
  expect_identical(sData$type, "S")
  expect_equal(pData$data$n, c(0, 1, 2))
  expect_equal(pData$data$rn, c(20, 5, 1))
  expect_equal(sData$data$n, c(1, 2, 3))
  expect_equal(sData$data$rn, c(8, 3, 1))
})

test_that("makePSData expands observations when counts are omitted", {
  pData = makePSData(n = c(0, 0, 0, 1, 2), type = "P")

  expect_s3_class(pData, "psData")
  expect_equal(pData$data$n, c(0, 1, 2))
  expect_equal(pData$data$rn, c(3, 1, 1))
})

test_that("readData reads a two-column CSV into psData", {
  csvPath = tempfile(fileext = ".csv")
  write.csv(
    data.frame(P = c(0, 1, 2), count = c(20, 5, 1)),
    csvPath,
    row.names = FALSE
  )

  pData = readData(csvPath)

  expect_s3_class(pData, "psData")
  expect_identical(pData$type, "P")
  expect_equal(pData$data$n, c(0, 1, 2))
  expect_equal(pData$data$rn, c(20, 5, 1))
})

test_that("readData rejects duplicate categories", {
  csvPath = tempfile(fileext = ".csv")
  write.csv(
    data.frame(S = c(1, 1), count = c(2, 3)),
    csvPath,
    row.names = FALSE
  )

  expect_error(readData(csvPath), "unique")
})
