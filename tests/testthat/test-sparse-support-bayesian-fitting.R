test_that("Pettard S data retain information for Bayesian fitting", {
  data(Ssurveys)
  pettard = Ssurveys$pettard

  expect_identical(pettard$type, "S")
  expect_equal(pettard$data$n, 1)
  expect_equal(pettard$data$rn, 6)
  expect_equal(modelObservationData(zetaModel(), pettard), 1)

  expect_true(is.finite(
    modelLogLikelihood(zetaModel(), list(shape = 2), pettard)
  ))
  expect_true(is.finite(
    modelLogLikelihood(logarithmicModel(), list(pi = 0.5), pettard)
  ))
  expect_true(is.finite(
    modelLogLikelihood(zizModel(), list(pi = 0.5, shape = 2), pettard)
  ))
})

test_that("Bayesian fitting proceeds for Pettard S data with a prior warning", {
  data(Ssurveys)
  pettard = Ssurveys$pettard

  zetaFit = NULL
  expect_warning(
    {
      zetaFit = fit(
        pettard,
        model = zetaModel(),
        method = "bayes",
        bayesOptions = list(posteriorMethod = "numerical")
      )
    },
    "only one occupied support value"
  )
  expect_s3_class(zetaFit, "psFit")
  expect_identical(zetaFit$method, "bayes")
  expect_true(is.finite(zetaFit$shape))

  logarithmicFit = NULL
  expect_warning(
    {
      logarithmicFit = fit(
        pettard,
        model = logarithmicModel(),
        method = "bayes",
        bayesOptions = list(posteriorMethod = "numerical")
      )
    },
    "only one occupied support value"
  )
  expect_s3_class(logarithmicFit, "psFit")
  expect_identical(logarithmicFit$method, "bayes")
  expect_true(is.finite(logarithmicFit$pi))
})

test_that("two-parameter Bayesian fitting can use Pettard S data", {
  skip_if_not_installed("cubature")
  data(Ssurveys)
  pettard = Ssurveys$pettard

  zizFit = NULL
  expect_warning(
    {
      zizFit = fit(
        pettard,
        model = zizModel(),
        method = "bayes",
        bayesOptions = list(posteriorMethod = "numerical"),
        tol = 1e-4,
        summaryGridSize = 21L
      )
    },
    "only one occupied support value"
  )

  expect_s3_class(zizFit, "psFit")
  expect_identical(zizFit$method, "bayes")
  expect_true(is.finite(zizFit$pi))
  expect_true(is.finite(zizFit$shape))
})

test_that("historical MLE singleton-support rejection is preserved", {
  data(Ssurveys)
  pettard = Ssurveys$pettard

  expect_error(
    fit(pettard, model = zetaModel(), method = "mle"),
    "at least one value higher than 1"
  )
  expect_error(
    fit(pettard, model = logarithmicModel(), method = "mle"),
    "at least one value higher than 1"
  )
  expect_error(
    fit(pettard, model = zizModel(), method = "mle"),
    "at least one value higher than 1"
  )
})
