test_that("posteriorInflation respects numerical posterior weights", {
  posterior = structure(
    list(
      method = "numerical",
      representation = list(
        pi = c(0.005, 0.02),
        marginalDensity = list(pi = c(40, 10)),
        dPi = 0.02
      ),
      model = "ziz"
    ),
    class = "psPosterior"
  )

  result = posteriorInflation(posterior, epsilon = 0.01)
  expect_equal(result$probBelow, 0.8)
  expect_equal(result$probAtOrAbove, 0.2)
})

test_that("posteriorInflation respects MCMC draws", {
  posterior = structure(
    list(
      method = "mcmc",
      representation = data.frame(pi = c(0.001, 0.005, 0.02, 0.03)),
      model = "ziz"
    ),
    class = "psPosterior"
  )

  result = posteriorInflation(posterior, epsilon = 0.01)
  expect_equal(result$probBelow, 0.5)
})

test_that("posteriorInflation respects importance weights", {
  posterior = structure(
    list(
      method = "importance",
      representation = list(
        samples = data.frame(
          pi = c(0.001, 0.02, 0.03),
          weight = c(0.6, 0.3, 0.1)
        )
      ),
      model = "ziz"
    ),
    class = "psPosterior"
  )

  result = posteriorInflation(posterior, epsilon = 0.01)
  expect_equal(result$probBelow, 0.6)
})

test_that("posteriorInflation uses the Laplace logit approximation", {
  posterior = structure(
    list(
      method = "laplace",
      representation = list(
        modeWorking = c(eta = -4, tau = 0),
        covarianceWorking = matrix(
          c(0.25, 0, 0, 0.5),
          nrow = 2,
          dimnames = list(c("eta", "tau"), c("eta", "tau"))
        )
      ),
      model = "ziz"
    ),
    class = "psPosterior"
  )

  expected = pnorm(qlogis(0.01), mean = -4, sd = 0.5)
  result = posteriorInflation(posterior, epsilon = 0.01)
  expect_equal(result$probBelow, expected)
})

test_that("posteriorInflation validates epsilon", {
  posterior = structure(
    list(
      method = "mcmc",
      representation = data.frame(pi = c(0.001, 0.02)),
      model = "ziz"
    ),
    class = "psPosterior"
  )

  expect_error(posteriorInflation(posterior, epsilon = 0), "strictly between")
  expect_error(posteriorInflation(posterior, epsilon = 1), "strictly between")
})

test_that("Roux summary reports the practical inflation probability", {
  data(Psurveys)
  fit = fitZIDist(
    Psurveys$roux,
    nterms = 4,
    method = "bayes",
    bayesOptions = list(posteriorMethod = "numerical"),
    nPiGrid = 41,
    nShapeGrid = 41
  )

  direct = posteriorInflation(fit, epsilon = 0.01)
  summarised = summary(fit, inflationEpsilon = 0.01)

  expect_s3_class(summarised, "summary.psPosterior")
  expect_equal(summarised$inflation, direct)
  expect_gte(direct$probBelow, 0)
  expect_lte(direct$probBelow, 1)
})
