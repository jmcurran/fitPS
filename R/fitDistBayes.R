#' Fit the plain zeta model with the MCMC Bayesian implementation.
#'
#' @param x An object of class `psData`.
#' @param prior A prior object created by `makePrior()`.
#' @param nterms Number of fitted P/S probability terms to retain.
#' @param ... MCMC controls passed to the posterior engine, including `shape0`,
#'   `nIter`, `nBurnIn`, and `silent`.
#' @return A Bayesian `psFit` object.
#' @keywords internal
#' @noRd
fitDistBayes = function(x, prior = makePrior(), nterms, ...) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData")
  }

  model = zetaModel()
  engine = mcmcPosteriorEngine()
  representation = fitPosterior(
    engine = engine,
    model = model,
    x = x,
    prior = prior,
    ...
  )

  chain = representation$value$chain
  densityEstimate = density(chain, from = prior$range[1])
  pdf = splinefun(densityEstimate$x, densityEstimate$y)
  varShape = unname(representation$value$variance["shape", "shape"])

  finaliseBayesianPsFit(
    model = model,
    engine = engine,
    representation = representation,
    x = x,
    nterms = nterms,
    fit = list(par = unname(representation$value$mean[["shape"]])),
    legacyFields = list(
      var.shape = varShape,
      chain = chain,
      pdf = pdf
    )
  )
}
