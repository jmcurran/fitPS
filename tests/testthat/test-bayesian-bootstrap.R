test_that("grouped Bayesian Bootstrap weights are positive and preserve total mass", {
  pData = makePSData(
    n = c(0, 1, 2, 3),
    count = c(40, 12, 4, 1),
    type = "P"
  )

  weights = drawGroupedBayesianBootstrapWeights(
    pData,
    B = 25,
    seed = 1103
  )

  expect_equal(dim(weights), c(25L, 4L))
  expect_true(all(weights > 0))
  expect_equal(
    rowSums(weights),
    rep(sum(pData$data$rn), 25),
    tolerance = 1e-12
  )
  expect_identical(colnames(weights), as.character(pData$data$n))
})

test_that("grouped Bayesian Bootstrap draws are reproducible under a seed", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(10, 5, 2),
    type = "P"
  )

  first = drawGroupedBayesianBootstrapWeights(pData, B = 20, seed = 1104)
  second = drawGroupedBayesianBootstrapWeights(pData, B = 20, seed = 1104)
  third = drawGroupedBayesianBootstrapWeights(pData, B = 20, seed = 1105)

  expect_identical(first, second)
  expect_false(identical(first, third))
})

test_that("grouped Dirichlet draws reproduce category-count mean weights", {
  pData = makePSData(
    n = c(0, 1, 2),
    count = c(12, 6, 2),
    type = "P"
  )

  weights = drawGroupedBayesianBootstrapWeights(
    pData,
    B = 4000,
    seed = 1106,
    sampleSizeScale = FALSE
  )

  expected = pData$data$rn / sum(pData$data$rn)
  names(expected) = as.character(pData$data$n)
  expect_equal(colMeans(weights), expected, tolerance = 0.015)
  expect_equal(rowSums(weights), rep(1, 4000), tolerance = 1e-12)
})

test_that("Bayesian Bootstrap weighted refits retain replicate dimensions and names", {
  pData = makePSData(
    n = c(0, 1, 2, 3),
    count = c(40, 12, 4, 1),
    type = "P"
  )

  result = bayesianBootstrapModel(
    x = pData,
    model = zetaModel(),
    B = 30,
    seed = 1107,
    nterms = 4,
    level = 0.9
  )

  expect_s3_class(result, "psBayesianBootstrap")
  expect_identical(result$method, "rubin")
  expect_identical(result$model, "zeta")
  expect_identical(result$seed, 1107L)
  expect_equal(nrow(result$replicates$parameters), 30)
  expect_equal(nrow(result$replicates$probabilities), 30)
  expect_identical(names(result$replicates$parameters), "shape")
  expect_identical(
    names(result$replicates$probabilities),
    c("P0", "P1", "P2", "P3")
  )
  expect_identical(result$parameters$parameter, "shape")
  expect_identical(result$probabilities$term, c("P0", "P1", "P2", "P3"))
  expect_equal(result$level, 0.9)
  expect_equal(
    result$diagnostics$nSuccessful + result$diagnostics$nFailed,
    30
  )
})

test_that("Bayesian Bootstrap model replicates are seeded reproducibly", {
  pData = makePSData(
    n = c(0, 1, 2, 3),
    count = c(30, 12, 5, 2),
    type = "P"
  )

  first = bayesianBootstrapModel(
    pData,
    model = logarithmicModel(),
    B = 20,
    seed = 1108,
    nterms = 3
  )
  second = bayesianBootstrapModel(
    pData,
    model = logarithmicModel(),
    B = 20,
    seed = 1108,
    nterms = 3
  )

  expect_identical(first$replicates, second$replicates)
  expect_identical(first$parameters, second$parameters)
  expect_identical(first$probabilities, second$probabilities)
})

test_that("Bayesian Bootstrap retains failed weighted fits as missing rows", {
  sData = makePSData(n = 1, count = 6, type = "S")

  result = bayesianBootstrapModel(
    x = sData,
    model = zetaModel(),
    B = 8,
    seed = 1109,
    nterms = 3
  )

  expect_equal(result$diagnostics$nRequested, 8)
  expect_equal(result$diagnostics$nSuccessful, 0)
  expect_equal(result$diagnostics$nFailed, 8)
  expect_equal(result$diagnostics$failureRate, 1)
  expect_equal(nrow(result$diagnostics$failures), 8)
  expect_true(all(is.na(result$replicates$parameters$shape)))
  expect_true(all(vapply(
    result$diagnostics$failures$message,
    function(message) grepl("at least one value higher than 1", message),
    logical(1L)
  )))
  expect_true(all(is.na(result$parameters$estimate)))
  expect_true(all(is.na(result$probabilities$estimate)))
})

test_that("Bayesian Bootstrap supports all built-in weighted MLE models", {
  pData = makePSData(
    n = c(0, 1, 2, 3),
    count = c(50, 15, 5, 2),
    type = "P"
  )
  models = list(
    zetaModel(),
    zizModel(),
    logarithmicModel()
  )

  for (model in models) {
    result = bayesianBootstrapModel(
      pData,
      model = model,
      B = 12,
      seed = 1110,
      nterms = 3
    )

    expect_s3_class(result, "psBayesianBootstrap")
    expect_equal(nrow(result$replicates$parameters), 12)
    expect_equal(nrow(result$replicates$probabilities), 12)
    expect_equal(
      result$diagnostics$nSuccessful + result$diagnostics$nFailed,
      12
    )
  }
})
