test_that("Roux data support harmonised probability workflows", {
  roux = Psurveys$roux

  mleFit = fitZIDist(roux, nterms = 4)
  expect_s3_class(mleFit, "psFit")
  expect_equal(length(fitted(mleFit, n = 4, type = "plugIn")), 4L)

  bayesFit = fitZIDist(
    roux,
    nterms = 4,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "numerical"),
    nPiGrid = 21,
    nShapeGrid = 21
  )

  expect_s3_class(bayesFit$posterior, "psPosterior")
  expect_equal(
    names(fitted(bayesFit, n = 4, type = "posteriorMean")),
    paste0("P", 0:3)
  )
  expect_true(all(posteriorProbs(bayesFit, n = 4)$lower >= 0))
  expect_true(all(posteriorProbs(bayesFit, n = 4)$upper <= 1))

  bootFit = fit(
    roux,
    model = zizModel(),
    method = "bootstrap",
    nterms = 4,
    B = 12,
    seed = 51511,
    silent = TRUE,
    parallel = FALSE
  )

  expect_s3_class(bootFit$bootstrap, "psBootstrap")
  expect_equal(
    names(fitted(bootFit, n = 4, type = "bootstrapMean")),
    paste0("P", 0:3)
  )
  expect_true(all(bootstrapProbs(bootFit, n = 4)$lower >= 0))
  expect_true(all(bootstrapProbs(bootFit, n = 4)$upper <= 1))
})

test_that("legacy plug-in probability behaviour remains available on Roux data", {
  fit = fitZIDist(Psurveys$roux, nterms = 5)

  expect_equal(fitted(fit), fitted(fit, type = "plugIn"))

  P = probfun(fit)
  expect_equal(unname(P(0:4)), unname(fitted(fit, n = 5, type = "plugIn")))
})
