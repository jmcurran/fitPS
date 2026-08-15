test_that("built-in models satisfy the generic MCMC model contract", {
  pData = makePSData(
    n = c(0, 1, 2, 3),
    count = c(20, 8, 3, 1),
    type = "P"
  )

  cases = list(
    list(
      model = zetaModel(),
      prior = makePrior(family = "uniform", range = c(1.1, 4)),
      proposalScale = 0.2,
      dots = list()
    ),
    list(
      model = logarithmicModel(),
      prior = makePrior(family = "uniform", range = c(0.001, 0.999)),
      proposalScale = 0.3,
      dots = list()
    ),
    list(
      model = zizModel(),
      prior = makePrior(family = "uniform", range = c(1.1, 4)),
      proposalScale = c(pi = 0.25, shape = 0.2),
      dots = list(shape1 = 2, shape2 = 5)
    )
  )

  for (case in cases) {
    args = c(
      list(
        model = case$model,
        engine = mcmcPosteriorEngine(),
        x = pData,
        prior = case$prior,
        nIter = 250L,
        nBurnIn = 100L,
        proposalScale = case$proposalScale,
        seed = 20260815
      ),
      case$dots
    )

    representation = do.call(fitMcmcPosteriorModel.psModel, args)

    expect_s3_class(representation, "mcmcPosteriorRepresentation")
    expect_true(isTRUE(representation$metadata$generic))
    expect_identical(
      names(representation$value$mean),
      modelParameterNames(case$model)
    )
    expect_true(all(is.finite(representation$value$mean)))
    expect_true(all(is.finite(representation$value$variance)))
    expect_true(is.finite(representation$metadata$acceptance))
    expect_gt(representation$metadata$acceptance, 0)
    expect_lt(representation$metadata$acceptance, 1)
  }
})

test_that("specialised built-in MCMC methods remain explicit compatibility paths", {
  expect_true(is.function(getS3method(
    "fitMcmcPosteriorModel",
    "zetaModel",
    optional = TRUE
  )))
  expect_true(is.function(getS3method(
    "fitMcmcPosteriorModel",
    "logarithmicModel",
    optional = TRUE
  )))
  expect_true(is.function(getS3method(
    "fitMcmcPosteriorModel",
    "zizModel",
    optional = TRUE
  )))
})
