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
  summary = summarisePosterior(engine, model, representation)
  pointEstimate = posteriorPointEstimate(engine, model, representation)

  shape = unname(pointEstimate["shape"])
  varShape = summary$sd[summary$parameter == "shape"]^2
  probabilityIndices = if (x$type == "P") {
    0:(nterms - 1L)
  } else {
    seq_len(nterms)
  }
  fittedMatrix = modelProbabilities(
    model = model,
    parameters = pointEstimate,
    n = probabilityIndices,
    type = x$type
  )
  fitted = as.numeric(fittedMatrix[1L, ])
  names(fitted) = colnames(fittedMatrix)

  chain = representation$value$chain
  densityEstimate = density(chain, from = prior$range[1])
  pdf = splinefun(densityEstimate$x, densityEstimate$y)

  result = list(
    psData = x,
    fit = list(par = shape),
    shape = shape,
    var.shape = varShape,
    fitted = fitted,
    chain = chain,
    pdf = pdf,
    posteriorRepresentation = representation,
    model = "zeta",
    method = "bayes"
  )

  class(result) = "psFit"
  result
}
