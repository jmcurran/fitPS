#' Posterior probability of practically negligible inflation
#'
#' Calculate the posterior probability that the zero/one-inflation parameter
#' `pi` is smaller than a user-specified practical threshold.
#'
#' @param object A Bayesian zero-inflated `psFit` object or its associated
#'   `psPosterior` object.
#' @param epsilon A number strictly between 0 and 1 defining the largest
#'   inflation probability regarded as practically negligible.
#' @param ... Additional arguments passed to methods.
#'
#' @return A data frame with the threshold and posterior probabilities below
#'   and at or above that threshold.
#'
#' @details
#' With the continuous beta prior used for `pi`, the posterior probability of
#' exactly `pi = 0` is zero. A more practically useful diagnostic is therefore
#' `Pr(pi < epsilon | data)`, where `epsilon` is chosen to represent an
#' inflation effect small enough to be negligible for the application.
#'
#' The calculation respects the stored posterior engine. Numerical fits
#' integrate the marginal posterior grid, MCMC fits use retained draws,
#' importance fits retain their importance weights, and Laplace fits evaluate
#' the implied Gaussian approximation on the logit scale.
#'
#' @examples
#' data(Psurveys)
#' roux = Psurveys$roux
#' if (interactive()) {
#'   fit = fitZIDist(
#'     roux,
#'     method = "bayes",
#'     bayesOptions = list(posteriorMethod = "numerical")
#'   )
#'   posteriorInflation(fit, epsilon = 0.01)
#' }
#'
#' @export
posteriorInflation = function(object, ...) {
  UseMethod("posteriorInflation")
}

#' @rdname posteriorInflation
#' @export
posteriorInflation.psFit = function(object, epsilon = 0.01, ...) {
  if (!identical(object$model, "ziz") ||
      !identical(object$method, "bayes") ||
      !inherits(object$posterior, "psPosterior")) {
    stop("posteriorInflation() is only available for Bayesian zero-inflated psFit objects")
  }

  posteriorInflation(object$posterior, epsilon = epsilon, ...)
}

#' @rdname posteriorInflation
#' @importFrom stats pnorm qlogis
#' @export
posteriorInflation.psPosterior = function(object, epsilon = 0.01, ...) {
  validateInflationEpsilon(epsilon)

  if (!is.null(object$model) && !identical(object$model, "ziz")) {
    stop("posteriorInflation() requires a zero-inflated posterior")
  }

  probabilityBelow = posteriorInflationProbability(
    object$representation,
    epsilon = epsilon
  )

  probabilityBelow = min(1, max(0, unname(probabilityBelow)))

  data.frame(
    epsilon = epsilon,
    probBelow = probabilityBelow,
    probAtOrAbove = 1 - probabilityBelow,
    posteriorMethod = object$method,
    stringsAsFactors = FALSE
  )
}

#' Validate the practical zero-inflation threshold used by posterior diagnostics.
#'
#' @param epsilon Practical threshold below which zero inflation is treated as negligible.
#' @return The validated threshold, invisibly.
#' @keywords internal
#' @noRd
validateInflationEpsilon = function(epsilon) {
  if (!is.numeric(epsilon) || length(epsilon) != 1L ||
      !is.finite(epsilon) || epsilon <= 0 || epsilon >= 1) {
    stop("epsilon must be one finite number strictly between 0 and 1")
  }

  invisible(epsilon)
}

#' Compute the posterior probability below the inflation threshold from a numerical grid.
#'
#' @param posteriorGrid Numerical posterior-grid representation.
#' @param epsilon Practical threshold below which zero inflation is treated as negligible.
#' @return A posterior probability between zero and one.
#' @keywords internal
#' @noRd
numericalInflationProbability = function(posteriorGrid, epsilon) {
  required = c("pi", "marginalDensity", "dPi")
  if (!is.list(posteriorGrid) || !all(required %in% names(posteriorGrid)) ||
      is.null(posteriorGrid$marginalDensity$pi)) {
    stop("numerical posterior representation is incomplete")
  }

  # With a continuous prior Pr(pi = 0 | data) is zero. The useful diagnostic is
  # therefore the posterior mass below a small practical threshold epsilon.
  selected = posteriorGrid$pi < epsilon
  sum(posteriorGrid$marginalDensity$pi[selected]) * posteriorGrid$dPi
}

#' Compute the posterior probability below the inflation threshold from MCMC draws.
#'
#' @param chain MCMC posterior draws.
#' @param epsilon Practical threshold below which zero inflation is treated as negligible.
#' @return A posterior probability between zero and one.
#' @keywords internal
#' @noRd
mcmcInflationProbability = function(chain, epsilon) {
  if (!is.data.frame(chain) || !"pi" %in% names(chain)) {
    stop("MCMC posterior representation does not contain pi draws")
  }

  # Estimate the practical no-inflation probability by the retained-draw
  # proportion below epsilon; exact zero has probability zero under the prior.
  mean(chain$pi < epsilon)
}

#' Compute the posterior probability below the inflation threshold from weighted importance samples.
#'
#' @param approximation Laplace or importance-sampling posterior representation.
#' @param epsilon Practical threshold below which zero inflation is treated as negligible.
#' @return A posterior probability between zero and one.
#' @keywords internal
#' @noRd
importanceInflationProbability = function(approximation, epsilon) {
  if (!is.list(approximation) || is.null(approximation$samples) ||
      !all(c("pi", "weight") %in% names(approximation$samples))) {
    stop("importance posterior representation does not contain weighted pi samples")
  }

  # Weighted posterior mass below epsilon estimates the practical no-inflation
  # probability for an importance-sampling representation.
  samples = approximation$samples
  sum(samples$weight[samples$pi < epsilon]) / sum(samples$weight)
}

#' Approximate the posterior probability below the inflation threshold from a Laplace approximation.
#'
#' @param approximation Laplace or importance-sampling posterior representation.
#' @param epsilon Practical threshold below which zero inflation is treated as negligible.
#' @return An approximate posterior probability between zero and one.
#' @keywords internal
#' @noRd
laplaceInflationProbability = function(approximation, epsilon) {
  if (!is.list(approximation) || is.null(approximation$modeWorking) ||
      is.null(approximation$covarianceWorking)) {
    stop("Laplace posterior representation is incomplete")
  }

  etaMean = unname(approximation$modeWorking[["eta"]])
  etaSd = sqrt(unname(approximation$covarianceWorking["eta", "eta"]))

  if (!is.finite(etaMean) || !is.finite(etaSd) || etaSd <= 0) {
    stop("Laplace approximation for pi is not valid")
  }

  # The Laplace approximation is Gaussian in eta = logit(pi), so transform the
  # natural-scale threshold to eta before evaluating its Gaussian CDF.
  pnorm(qlogis(epsilon), mean = etaMean, sd = etaSd)
}

#' Compute negligible-inflation probability through posterior representation dispatch.
#'
#' @param representation An internal `psPosteriorRepresentation` object.
#' @param epsilon Practical inflation threshold.
#' @param ... Additional representation-specific controls.
#' @return Posterior probability that `pi < epsilon`.
#' @keywords internal
#' @noRd
posteriorInflationProbability = function(representation, epsilon, ...) {
  UseMethod("posteriorInflationProbability")
}

#' @rdname posteriorInflationProbability
#' @keywords internal
#' @exportS3Method posteriorInflationProbability psPosteriorRepresentation
#' @noRd
posteriorInflationProbability.psPosteriorRepresentation = function(representation,
                                                                    epsilon,
                                                                    ...) {
  stop(
    "Inflation probability is not implemented for posterior representation '",
    class(representation)[1L],
    "'",
    call. = FALSE
  )
}

#' @rdname posteriorInflationProbability
#' @keywords internal
#' @exportS3Method posteriorInflationProbability numericalPosteriorRepresentation
#' @noRd
posteriorInflationProbability.numericalPosteriorRepresentation = function(representation,
                                                                           epsilon,
                                                                           ...) {
  validateInflationEpsilon(epsilon)
  posteriorGrid = representation$value$posteriorGrid
  numericalInflationProbability(posteriorGrid, epsilon)
}

#' @rdname posteriorInflationProbability
#' @keywords internal
#' @exportS3Method posteriorInflationProbability mcmcPosteriorRepresentation
#' @noRd
posteriorInflationProbability.mcmcPosteriorRepresentation = function(representation,
                                                                      epsilon,
                                                                      ...) {
  validateInflationEpsilon(epsilon)
  mcmcInflationProbability(representation$value$chain, epsilon)
}

#' @rdname posteriorInflationProbability
#' @keywords internal
#' @exportS3Method posteriorInflationProbability importancePosteriorRepresentation
#' @noRd
posteriorInflationProbability.importancePosteriorRepresentation = function(representation,
                                                                            epsilon,
                                                                            ...) {
  validateInflationEpsilon(epsilon)
  importanceInflationProbability(representation$value$approximation, epsilon)
}

#' @rdname posteriorInflationProbability
#' @keywords internal
#' @exportS3Method posteriorInflationProbability laplacePosteriorRepresentation
#' @noRd
posteriorInflationProbability.laplacePosteriorRepresentation = function(representation,
                                                                         epsilon,
                                                                         ...) {
  validateInflationEpsilon(epsilon)
  laplaceInflationProbability(representation$value$approximation, epsilon)
}
