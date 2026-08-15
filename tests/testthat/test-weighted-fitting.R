test_that("weighted fitting preserves ordinary built-in MLEs at observed counts", {
  pData = makePSData(
    n = c(0, 1, 2, 3),
    count = c(40, 12, 4, 1),
    type = "P"
  )
  models = list(
    zetaModel(),
    zizModel(),
    logarithmicModel()
  )

  for (model in models) {
    ordinaryFit = fit(pData, model = model, nterms = 5)
    weightedFit = fitWeightedModel(
      pData,
      model = model,
      weights = pData$data$rn,
      nterms = 5
    )

    ordinaryParameters = fitModelParameters(ordinaryFit, model)
    expect_equal(
      unname(weightedFit$parameters),
      unname(ordinaryParameters),
      tolerance = 1e-6
    )
    expect_equal(weightedFit$fitted, ordinaryFit$fitted, tolerance = 1e-6)
    expect_identical(weightedFit$psData, pData)
  }
})

test_that("weighted fitting accepts positive fractional category weights", {
  pData = makePSData(
    n = c(0, 1, 2, 3),
    count = c(40, 12, 4, 1),
    type = "P"
  )
  weights = c(38.25, 13.5, 4.75, 0.5)
  models = list(
    zetaModel(),
    zizModel(),
    logarithmicModel()
  )

  for (model in models) {
    weightedFit = fitWeightedModel(
      pData,
      model = model,
      weights = weights,
      nterms = 5
    )

    expect_s3_class(weightedFit, "psWeightedMle")
    expect_identical(weightedFit$weights, weights)
    expect_true(all(is.finite(weightedFit$parameters)))
    expect_true(all(is.finite(weightedFit$fitted)))
    expect_equal(sum(weightedFit$psData$data$rn), sum(pData$data$rn))
  }
})

test_that("weighted likelihood is invariant to common positive weight scaling", {
  pData = makePSData(
    n = c(0, 1, 2, 3),
    count = c(40, 12, 4, 1),
    type = "P"
  )
  weights = c(0.61, 0.23, 0.11, 0.05)
  models = list(
    zetaModel(),
    zizModel(),
    logarithmicModel()
  )

  for (model in models) {
    unitScale = fitWeightedModel(pData, model, weights = weights, nterms = 4)
    sampleScale = fitWeightedModel(
      pData,
      model,
      weights = sum(pData$data$rn) * weights,
      nterms = 4
    )

    expect_equal(
      unitScale$parameters,
      sampleScale$parameters,
      tolerance = 1e-5
    )
    expect_equal(unitScale$fitted, sampleScale$fitted, tolerance = 1e-5)
  }
})

test_that("weighted fitting rejects invalid weights without changing psData semantics", {
  pData = makePSData(n = c(0, 1, 2), count = c(8, 3, 1), type = "P")

  expect_error(
    fitWeightedModel(pData, zetaModel(), weights = c(1, 2)),
    "aligned with the occupied support rows"
  )
  expect_error(
    fitWeightedModel(pData, zetaModel(), weights = c(1, 0, 2)),
    "finite positive values"
  )
  expect_error(
    fitWeightedModel(pData, zetaModel(), weights = c(1, NA, 2)),
    "finite positive values"
  )

  expect_warning(
    makePSData(
      n = c(0, 1, 2),
      count = c(1.2, 2.2, 3.2),
      type = "P"
    ),
    "Non-integer input detected"
  )
  rounded = suppressWarnings(
    makePSData(
      n = c(0, 1, 2),
      count = c(1.2, 2.2, 3.2),
      type = "P"
    )
  )
  expect_identical(rounded$data$rn, c(1, 2, 3))
})

test_that("weighted fitting retains the historical sparse-support MLE rule", {
  sData = makePSData(n = 1, count = 6, type = "S")

  expect_error(
    fitWeightedModel(sData, zetaModel(), weights = 6),
    "at least one value higher than 1"
  )
})
