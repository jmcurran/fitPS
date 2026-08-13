#' Fit a numerical posterior for the zero-inflated zeta model.
#'
#' The ZIZ model retains its existing rectangular two-dimensional grid
#' calculation. The numerical engine wraps that representation in the common
#' Stage 6 posterior contract rather than duplicating the integration logic.
#'
#' @param model A `zizModel` descriptor.
#' @param engine A `numericalPosteriorEngine` object.
#' @param x An object of class `psData`.
#' @param prior A prior object created by `makePrior()`.
#' @param shape1 First beta-prior shape parameter for `pi`.
#' @param shape2 Second beta-prior shape parameter for `pi`.
#' @param nPiGrid Number of grid points for `pi`.
#' @param nShapeGrid Number of grid points for `shape`.
#' @param ... Additional numerical controls reserved for future use.
#' @return A `numericalPosteriorRepresentation` object.
#' @keywords internal
#' @noRd
fitNumericalPosteriorModel.zizModel = function(model,
                                                engine,
                                                x,
                                                prior,
                                                shape1 = 1,
                                                shape2 = 1,
                                                nPiGrid = 101,
                                                nShapeGrid = 101,
                                                ...) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData")
  }
  if (!inherits(prior, "psPrior")) {
    stop("prior must be an object of class psPrior")
  }

  modelObservationData(model, x)
  posteriorGrid = makeZizPosteriorGrid(
    x = x,
    prior = prior,
    shape1 = shape1,
    shape2 = shape2,
    nPiGrid = nPiGrid,
    nShapeGrid = nShapeGrid
  )

  newPsPosteriorRepresentation(
    engine = engine,
    value = list(
      posteriorGrid = posteriorGrid,
      mean = posteriorGrid$mean,
      variance = posteriorGrid$varCov
    ),
    metadata = list(
      model = model$model,
      nPiGrid = length(posteriorGrid$pi),
      nShapeGrid = length(posteriorGrid$shape),
      normalizingConstant = posteriorGrid$normalizingConstant
    )
  )
}
