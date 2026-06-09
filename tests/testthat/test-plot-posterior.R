test_that("plotPosterior returns density data for zeta MCMC samples", {
  fit = list(
    shape = 2.4,
    fit = list(par = 2.4),
    chain = seq(1.8, 3.0, length.out = 200),
    model = "zeta",
    method = "bayes"
  )
  class(fit) = "psFit"

  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)

  result = plotPosterior(fit)

  expect_s3_class(result, "data.frame")
  expect_named(result, c("x", "density"))
  expect_equal(attr(result, "estimate"), 2.4)
  expect_length(attr(result, "interval"), 2)
})

test_that("plotPosterior supports zero-inflated pi and shape samples", {
  chain = data.frame(
    pi = seq(0.2, 0.7, length.out = 200),
    shape = seq(1.5, 3.5, length.out = 200)
  )
  fit = list(
    pi = 0.45,
    shape = 2.5,
    fit = list(par = c(pi = 0.45, shape = 2.5)),
    chain = chain,
    model = "ziz",
    method = "bayes"
  )
  class(fit) = "psFit"

  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)

  piResult = plotPosterior(fit, parameter = "pi")
  shapeResult = plotPosterior(fit, parameter = "shape")

  expect_s3_class(piResult, "data.frame")
  expect_s3_class(shapeResult, "data.frame")
  expect_equal(attr(piResult, "estimate"), 0.45)
  expect_equal(attr(shapeResult, "estimate"), 2.5)
})

test_that("plotPosterior supports stored numerical posterior density", {
  fit = list(
    shape = 2.1,
    fit = list(par = 2.1),
    var.shape = 0.04,
    pdf = function(x) {
      exp(-0.5 * ((x - 2.1) / 0.2)^2)
    },
    model = "zeta",
    method = "integrate"
  )
  class(fit) = "psFit"

  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)

  result = plotPosterior(fit, nGrid = 128)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 128)
  expect_equal(attr(result, "estimate"), 2.1)
  expect_true(all(attr(result, "interval") > 1))
})

test_that("plotPosterior rejects non-Bayesian fits and unavailable parameters", {
  fit = list(
    shape = 2.4,
    model = "zeta",
    method = "mle"
  )
  class(fit) = "psFit"

  expect_error(plotPosterior(fit), "Bayesian")

  bayesFit = fit
  bayesFit$method = "bayes"
  expect_error(plotPosterior(bayesFit, parameter = "pi"), "zero-inflated")
})
