#' Fit the plain zeta model using one-dimensional numerical posterior integration.
#'
#' @param x An object of class `psData`.
#' @param prior A prior object created by `makePrior()`.
#' @param nterms Number of fitted P/S probability terms to retain.
#' @param ... Additional arguments passed to the numerical posterior engine.
#' @return A Bayesian `psFit` object.
#' @keywords internal
#' @noRd
fitDistBayesIntegrate = function(x, prior = makePrior(), nterms, ...) {
  model = zetaModel()
  engine = numericalPosteriorEngine()
  representation = fitPosterior(
    engine = engine,
    model = model,
    x = x,
    prior = prior,
    ...
  )
  posteriorSummary = summarisePosterior(
    engine = engine,
    model = model,
    representation = representation
  )
  pointEstimate = posteriorPointEstimate(
    engine = engine,
    model = model,
    representation = representation
  )

  shape = unname(pointEstimate[["shape"]])
  shapeSd = posteriorSummary$sd[posteriorSummary$parameter == "shape"]
  varShape = unname(shapeSd^2)

  probabilityIndices = if (x$type == "P") {
    0:(nterms - 1L)
  } else {
    seq_len(nterms)
  }
  fittedMatrix = modelProbabilities(
    model = model,
    parameters = list(shape = shape),
    n = probabilityIndices,
    type = x$type
  )
  fitted = as.numeric(fittedMatrix[1L, ])
  names(fitted) = colnames(fittedMatrix)

  result = list(
    psData = x,
    fit = list(),
    shape = shape,
    var.shape = varShape,
    fitted = fitted,
    pdf = representation$value$density,
    posteriorRepresentation = representation,
    model = "zeta",
    method = "integrate"
  )
  class(result) = "psFit"
  result
}
