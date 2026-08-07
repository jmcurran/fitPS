test_that("zetaProbabilities maps P and S terms on the standard shape scale", {
  shape = c(2, 2.5)

  pActual = zetaProbabilities(shape, n = 0:2, type = "P")
  pExpected = cbind(
    P0 = dzetaStandard(1, shape = shape),
    P1 = dzetaStandard(2, shape = shape),
    P2 = dzetaStandard(3, shape = shape)
  )
  expect_equal(pActual, pExpected)

  sActual = zetaProbabilities(shape, n = 1:3, type = "S")
  sExpected = cbind(
    S1 = dzetaStandard(1, shape = shape),
    S2 = dzetaStandard(2, shape = shape),
    S3 = dzetaStandard(3, shape = shape)
  )
  expect_equal(sActual, sExpected)
})

test_that("zetaProbabilities validates term indices", {
  expect_error(zetaProbabilities(2, -1, "P"), "non-negative")
  expect_error(zetaProbabilities(2, 0, "S"), "positive")
  expect_error(zetaProbabilities(1, 0, "P"), "greater than 1")
})
