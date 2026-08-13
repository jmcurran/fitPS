#' Fit the plain zeta model using one-dimensional numerical posterior integration.
#'
#' @param x An object of class `psData`.
#' @param prior A prior object created by `makePrior()`.
#' @param nterms Number of fitted P/S probability terms to retain.
#' @param level Equal-tailed credible interval level used for the common
#'   `psPosterior` probability summaries.
#' @param ... Additional arguments passed to the numerical posterior engine.
#' @return A Bayesian `psFit` object.
#' @keywords internal
#' @noRd
fitDistBayesIntegrate = function(x,
                                  prior = makePrior(),
                                  nterms,
                                  level = 0.95,
                                  ...) {
  model = zetaModel()
  engine = numericalPosteriorEngine()
  representation = fitPosterior(
    engine = engine,
    model = model,
    x = x,
    prior = prior,
    ...
  )

  varShape = unname(representation$value$variance["shape", "shape"])

  finaliseBayesianPsFit(
    model = model,
    engine = engine,
    representation = representation,
    x = x,
    nterms = nterms,
    level = level,
    fit = list(),
    legacyFields = list(
      var.shape = varShape,
      pdf = representation$value$density
    )
  )
}
