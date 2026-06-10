test_that("normaliseBayesMethod translates legacy Bayesian aliases", {
  expect_equal(
    normaliseBayesMethod("bayes")$method,
    "bayes"
  )
  expect_equal(
    normaliseBayesMethod("mle")$method,
    "mle"
  )

  expect_warning(
    info <- normaliseBayesMethod("integrate"),
    "deprecated"
  )
  expect_equal(info$method, "bayes")
  expect_equal(info$bayesOptions$posteriorMethod, "numerical")

  expect_warning(
    info <- normaliseBayesMethod("mcmc"),
    "deprecated"
  )
  expect_equal(info$method, "bayes")
  expect_equal(info$bayesOptions$posteriorMethod, "mcmc")
})

test_that("normaliseBayesMethod rejects conflicting legacy options", {
  expect_error(
    normaliseBayesMethod(
      "integrate",
      bayesOptions = list(posteriorMethod = "mcmc")
    ),
    "conflicts"
  )
})

test_that("Bayesian API defaults to numerical posterior method", {
  options = normaliseBayesOptions()

  expect_equal(options$posteriorMethod, "numerical")
  expect_true(is.list(options$prior))
})
