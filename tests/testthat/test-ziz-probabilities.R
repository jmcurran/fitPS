test_that("zizProbabilities maps P terms using the standard shape scale", {
  pi = c(0.2, 0.4)
  shape = c(2, 3)
  n = 0:2

  actual = zizProbabilities(pi, shape, n, type = "P")

  expected = cbind(
    P0 = pi + (1 - pi) * dzetaStandard(1, shape = shape),
    P1 = (1 - pi) * dzetaStandard(2, shape = shape),
    P2 = (1 - pi) * dzetaStandard(3, shape = shape)
  )

  expect_equal(actual, expected)
})

test_that("zizProbabilities maps S terms and names them correctly", {
  pi = c(0.15, 0.35)
  shape = 2.5
  n = 1:3

  actual = zizProbabilities(pi, shape, n, type = "S")

  expected = cbind(
    S1 = pi + (1 - pi) * dzetaStandard(1, shape = shape),
    S2 = (1 - pi) * dzetaStandard(2, shape = shape),
    S3 = (1 - pi) * dzetaStandard(3, shape = shape)
  )

  expect_equal(actual, expected)
  expect_equal(colnames(actual), paste0("S", n))
})

test_that("zizProbabilities supports scalar parameter recycling", {
  actual = zizProbabilities(
    pi = 0.25,
    shape = c(2, 2.5, 3),
    n = c(0, 1),
    type = "P"
  )

  expect_equal(dim(actual), c(3L, 2L))
  expect_equal(colnames(actual), c("P0", "P1"))
})

test_that("zizProbabilities approaches a complete distribution", {
  pProbabilities = zizProbabilities(
    pi = 0.3,
    shape = 2.5,
    n = 0:9999,
    type = "P"
  )
  sProbabilities = zizProbabilities(
    pi = 0.3,
    shape = 2.5,
    n = 1:10000,
    type = "S"
  )

  expect_equal(sum(pProbabilities), 1, tolerance = 1e-6)
  expect_equal(sum(sProbabilities), 1, tolerance = 1e-6)
})

test_that("zizProbabilities validates parameters and term indices", {
  expect_error(
    zizProbabilities(-0.1, 2, 0, "P"),
    "pi must contain"
  )
  expect_error(
    zizProbabilities(0.2, 1, 0, "P"),
    "greater than 1"
  )
  expect_error(
    zizProbabilities(c(0.2, 0.3), c(2, 2.5, 3), 0, "P"),
    "equal lengths"
  )
  expect_error(
    zizProbabilities(0.2, 2, -1, "P"),
    "non-negative"
  )
  expect_error(
    zizProbabilities(0.2, 2, 0, "S"),
    "positive"
  )
})
