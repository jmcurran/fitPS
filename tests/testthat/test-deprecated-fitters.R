test_that("legacy fitters warn and direct users to fit", {
  oldOption = getOption("fitPS.deprecationWarnings")
  on.exit(options(fitPS.deprecationWarnings = oldOption), add = TRUE)
  options(fitPS.deprecationWarnings = TRUE)

  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")

  expect_warning(
    fitDist(pData, nterms = 4),
    regexp = "fitDist\\(\\) is deprecated; use fit\\(x, model = zetaModel\\(\\)\\) instead"
  )
  expect_warning(
    fitZIDist(pData, nterms = 4),
    regexp = "fitZIDist\\(\\) is deprecated; use fit\\(x, model = zizModel\\(\\)\\) instead"
  )
  expect_warning(
    fitlogDist(pData, nterms = 4),
    regexp = "fitlogDist\\(\\) is deprecated; use fit\\(x, model = logarithmicModel\\(\\)\\) instead"
  )
})


test_that("legacy fitters retain model metadata and delegate through fit", {
  oldOption = getOption("fitPS.deprecationWarnings")
  on.exit(options(fitPS.deprecationWarnings = oldOption), add = TRUE)
  options(fitPS.deprecationWarnings = FALSE)

  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")

  legacyZeta = fitDist(pData, nterms = 4, start = 2)
  genericZeta = fit(pData, model = zetaModel(), nterms = 4, start = 2)
  expect_identical(legacyZeta$model, "zeta")
  expect_equal(legacyZeta$shape, genericZeta$shape)
  expect_equal(fitted(legacyZeta), fitted(genericZeta))

  legacyZiz = fitZIDist(pData, nterms = 4, start = c(0.25, 2))
  genericZiz = fit(pData, model = zizModel(), nterms = 4, start = c(0.25, 2))
  expect_identical(legacyZiz$model, "ziz")
  expect_equal(c(legacyZiz$pi, legacyZiz$shape), c(genericZiz$pi, genericZiz$shape))
  expect_equal(fitted(legacyZiz), fitted(genericZiz))

  legacyLogarithmic = fitlogDist(pData, nterms = 4, start = 0.5)
  genericLogarithmic = fit(
    pData,
    model = logarithmicModel(),
    nterms = 4,
    start = 0.5
  )
  expect_identical(legacyLogarithmic$model, "logarithmic")
  expect_equal(legacyLogarithmic$pi, genericLogarithmic$pi)
  expect_equal(fitted(legacyLogarithmic), fitted(genericLogarithmic))
})


test_that("bootstrapFit warns and directs users to fit bootstrap inference", {
  oldOption = getOption("fitPS.deprecationWarnings")
  on.exit(options(fitPS.deprecationWarnings = oldOption), add = TRUE)
  options(fitPS.deprecationWarnings = TRUE)

  pData = makePSData(n = c(0, 1, 2), count = c(20, 5, 1), type = "P")
  mleFit = fit(pData, model = zizModel(), nterms = 3)

  expect_warning(
    bootstrapFit(
      mleFit,
      B = 5,
      seed = 1115,
      silent = TRUE,
      parallel = FALSE
    ),
    regexp = paste0(
      "bootstrapFit\\(\\) is deprecated; use ",
      "fit\\(x, model = \\.\\.\\., method = \\\"bootstrap\\\"\\) instead"
    )
  )
})
