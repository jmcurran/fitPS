test_that("posterior grid CDF uses cumulative Simpson integration", {
  grid = seq(0, 1, length.out = 101)
  values = rep(1, length(grid))

  functions = makePosteriorGridFunctions(grid, values)

  expect_identical(functions$integrationRule, "simpson")
  expect_equal(functions$area, 1, tolerance = 1e-12)
  expect_equal(functions$cdf(c(0, 0.25, 0.5, 1)), c(0, 0.25, 0.5, 1), tolerance = 1e-12)
  expect_equal(functions$cumulative$cumulative[[1L]], 0)
  expect_equal(tail(functions$cumulative$cumulative, 1L), 1)
  expect_true(all(diff(functions$cumulative$cumulative) >= 0))

  quadratic = cumulativePosteriorIntegral(grid, grid^2)
  expect_identical(quadratic$method, "simpson")
  expect_equal(quadratic$area, 1 / 3, tolerance = 1e-12)
  expect_equal(quadratic$cumulative[seq(1, 101, by = 10)], grid[seq(1, 101, by = 10)]^3 / 3, tolerance = 1e-12)
})

test_that("numerical posterior contours use cumulative stored grid mass", {
  x = rep(c(0.2, 0.5, 0.8), each = 3)
  y = rep(c(1.5, 2.0, 2.5), times = 3)
  weights = c(
    0.01, 0.04, 0.01,
    0.04, 0.70, 0.08,
    0.01, 0.08, 0.03
  )
  weights = weights / sum(weights)

  region = makeGridUncertaintyRegion(
    grid = list(
      parameters = data.frame(pi = x, shape = y),
      weights = weights
    ),
    parameters = c("pi", "shape"),
    level = c(0.80, 0.95)
  )

  expect_equal(tail(region$cumulativeMass$cumulative, 1L), 1)
  expect_true(all(diff(region$cumulativeMass$cumulative) >= 0))
  expect_equal(names(region$thresholds), NULL)
  expect_length(region$thresholds, 2L)
  expect_gte(region$thresholds[[1L]], region$thresholds[[2L]])
})

test_that("one-dimensional replicate uncertainty returns density and intervals", {
  replicates = data.frame(shape = seq(1.2, 3.2, length.out = 200))
  parameterSummary = data.frame(
    parameter = "shape",
    estimate = mean(replicates$shape),
    sd = sd(replicates$shape),
    lower = quantile(replicates$shape, 0.025, names = FALSE),
    upper = quantile(replicates$shape, 0.975, names = FALSE),
    level = 0.95
  )
  probabilitySummary = data.frame(
    term = "P0",
    estimate = 0.5,
    sd = 0.1,
    lower = 0.3,
    upper = 0.7,
    level = 0.95,
    bootstrapMethod = "nonparametric"
  )
  bootstrap = newPsBootstrap(
    method = "nonparametric",
    parameters = parameterSummary,
    probabilities = probabilitySummary,
    replicates = replicates,
    level = 0.95
  )

  plotFile = tempfile(fileext = ".pdf")
  grDevices::pdf(plotFile)
  on.exit(grDevices::dev.off(), add = TRUE)

  result = plotUncertainty(bootstrap, level = c(0.80, 0.95))

  expect_identical(result$dimension, 1L)
  expect_identical(result$parameters, "shape")
  expect_equal(dim(result$intervals), c(2L, 2L))
  expect_true(all(c("x", "density") %in% names(result$density)))
})

test_that("Bayesian Bootstrap uncertainty uses stored parameter replicates", {
  replicates = expand.grid(
    pi = seq(0.05, 0.45, length.out = 15),
    shape = seq(1.4, 2.4, length.out = 15)
  )
  parameterSummary = data.frame(
    parameter = c("pi", "shape"),
    estimate = c(mean(replicates$pi), mean(replicates$shape)),
    sd = c(sd(replicates$pi), sd(replicates$shape)),
    lower = c(0.06, 1.50),
    upper = c(0.44, 1.75),
    level = c(0.95, 0.95)
  )
  probabilitySummary = data.frame(
    term = "P0",
    estimate = 0.5,
    sd = 0.1,
    lower = 0.3,
    upper = 0.7,
    level = 0.95
  )
  object = newPsBayesianBootstrap(
    method = "rubin",
    parameters = parameterSummary,
    probabilities = probabilitySummary,
    replicates = list(
      parameters = replicates,
      probabilities = data.frame(P0 = seq(0.3, 0.7, length.out = nrow(replicates)))
    ),
    level = 0.95,
    diagnostics = list(
      nRequested = nrow(replicates),
      nSuccessful = nrow(replicates),
      nFailed = 0L
    ),
    seed = 123,
    model = "ziz"
  )

  plotFile = tempfile(fileext = ".pdf")
  grDevices::pdf(plotFile)
  on.exit(grDevices::dev.off(), add = TRUE)

  result = plotUncertainty(
    object,
    level = c(0.80, 0.95),
    showPoints = FALSE
  )

  expect_identical(result$dimension, 2L)
  expect_identical(result$parameters, c("pi", "shape"))
  expect_true(length(result$contours) >= 1L)
  expect_true(all(vapply(result$contours, function(x) {
    x$level %in% c(0.80, 0.95)
  }, logical(1L))))
})

test_that("plotUncertainty validates higher-dimensional parameter selection", {
  expect_error(
    resolveUncertaintyParameters(c("a", "b", "c")),
    "more than two parameters"
  )
  expect_identical(
    resolveUncertaintyParameters(c("a", "b", "c"), c("a", "c")),
    c("a", "c")
  )
})

test_that("Laplace uncertainty contours use Gaussian probability radii", {
  contours = makeGaussianUncertaintyContours(
    centre = c(pi = 0.2, shape = 2),
    covariance = diag(c(0.01, 0.04)),
    level = c(0.80, 0.95),
    parameters = c("pi", "shape")
  )

  expect_length(contours, 2L)
  expect_equal(vapply(contours, `[[`, numeric(1L), "level"), c(0.80, 0.95))
  expect_true(all(vapply(contours, function(contour) {
    length(contour$x) == 241L && length(contour$y) == 241L
  }, logical(1L))))
})

test_that("weighted sample uncertainty accepts stored importance weights", {
  values = expand.grid(
    pi = seq(0.05, 0.45, length.out = 12),
    shape = seq(1.4, 2.4, length.out = 12)
  )
  weights = seq_len(nrow(values))

  region = makeSampleUncertaintyRegion(
    values = values,
    parameters = c("pi", "shape"),
    level = c(0.80, 0.95),
    weights = weights
  )

  expect_true(length(region$contours) >= 1L)
})
