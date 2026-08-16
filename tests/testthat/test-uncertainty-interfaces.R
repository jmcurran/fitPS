test_that("legacy bootstrap uncertainty interface directs users to fit and plotUncertainty", {
  oldOption = getOption("fitPS.deprecationWarnings")
  on.exit(options(fitPS.deprecationWarnings = oldOption), add = TRUE)
  options(fitPS.deprecationWarnings = TRUE)

  unsupported = structure(list(), class = "unsupportedUncertaintyInput")

  expect_warning(
    try(bootCI(unsupported), silent = TRUE),
    "bootCI\\(\\) is deprecated"
  )
})

test_that("legacy uncertainty deprecation warnings can be disabled", {
  oldOption = getOption("fitPS.deprecationWarnings")
  on.exit(options(fitPS.deprecationWarnings = oldOption), add = TRUE)
  options(fitPS.deprecationWarnings = FALSE)

  expect_silent(signalDeprecatedInterface("bootCI", "plotUncertainty()"))
  expect_silent(signalDeprecatedFitter("bootstrapFit", "fit(..., method = \"bootstrap\")"))
})

test_that("credint remains a Bayesian uncertainty extractor", {
  parameterReplicates = data.frame(
    pi = seq(0.1, 0.3, length.out = 60),
    shape = seq(1.8, 2.8, length.out = 60) + 0.08 * sin(seq(0, 6 * pi, length.out = 60))
  )
  probabilityReplicates = data.frame(P0 = rep(0.2, 60))
  parameterSummary = data.frame(
    parameter = c("pi", "shape"),
    estimate = c(0.2, 2.3),
    sd = c(0.05, 0.2),
    lower = c(0.1, 1.8),
    upper = c(0.3, 2.8),
    level = c(0.95, 0.95)
  )
  probabilitySummary = data.frame(
    term = "P0",
    estimate = 0.2,
    sd = 0.01,
    lower = 0.18,
    upper = 0.22,
    level = 0.95
  )
  object = newPsBayesianBootstrap(
    method = "Rubin Bayesian Bootstrap",
    parameters = parameterSummary,
    probabilities = probabilitySummary,
    replicates = list(
      parameters = parameterReplicates,
      probabilities = probabilityReplicates
    ),
    level = 0.95,
    diagnostics = list(),
    model = "ziz"
  )

  interval = credint(object, level = 0.95, parameters = "pi")
  expect_named(interval, c("lower", "upper"))
  expect_true(interval[[1L]] < interval[[2L]])

  region = credint(object, level = c(0.80, 0.95), silent = TRUE)
  expect_named(region, c("80%", "95%"))
  expect_true(all(vapply(region, length, integer(1L)) >= 1L))
})
