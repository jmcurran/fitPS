#' Draw from the Gaussian importance proposal in ZIZ working coordinates.
#'
#' @param mean Finite length-two proposal mean for `(eta, tau)`.
#' @param covariance Positive-definite 2 by 2 proposal covariance matrix.
#' @param n Number of proposal draws.
#' @return A numeric matrix with columns `eta` and `tau`.
#' @keywords internal
#' @noRd
#' @importFrom stats rnorm
makeZizProposalDraws = function(mean, covariance, n) {
  mean = unname(mean)

  if (!is.numeric(mean) || length(mean) != 2L || any(!is.finite(mean))) {
    stop("proposal mean must be a finite numeric vector of length two")
  }

  if (!is.matrix(covariance) || !all(dim(covariance) == c(2L, 2L))) {
    stop("proposal covariance must be a 2 by 2 matrix")
  }

  if (any(!is.finite(covariance))) {
    stop("proposal covariance must be finite")
  }

  covariance = (covariance + t(covariance)) / 2
  cholCovariance = tryCatch(
    chol(covariance),
    error = function(e) {
      NULL
    }
  )

  if (is.null(cholCovariance)) {
    stop("proposal covariance must be positive definite")
  }

  standardDraws = matrix(rnorm(n * 2L), ncol = 2L)
  draws = sweep(standardDraws %*% cholCovariance, 2L, mean, `+`)
  colnames(draws) = c("eta", "tau")
  draws
}

#' Evaluate the bivariate normal importance-proposal log density in working coordinates.
#'
#' @param working Numeric vector or matrix of working parameters `(eta, tau)`.
#' @param mean Mean vector used by the Gaussian proposal or weighted covariance calculation.
#' @param covariance Covariance matrix used by the Gaussian proposal.
#' @return Numeric log proposal densities.
#' @keywords internal
#' @noRd
zizProposalLogDensity = function(working, mean, covariance) {
  working = matrix(working, ncol = 2L)
  mean = unname(mean)
  covariance = (covariance + t(covariance)) / 2
  cholCovariance = chol(covariance)
  centered = sweep(working, 2L, mean, `-`)
  solved = backsolve(cholCovariance, t(centered), transpose = TRUE)
  quadratic = colSums(solved^2)
  logDeterminant = 2 * sum(log(diag(cholCovariance)))

  -log(2 * pi) - 0.5 * logDeterminant - 0.5 * quadratic
}

#' Compute a probability-weighted covariance matrix.
#'
#' @param values Matrix of values whose weighted covariance is required.
#' @param weights Non-negative normalised or normalisable weights.
#' @param mean Mean vector used by the Gaussian proposal or weighted covariance calculation.
#' @return A weighted covariance matrix.
#' @keywords internal
#' @noRd
weightedCovariance = function(values, weights, mean) {
  centered = sweep(values, 2L, mean, `-`)
  t(centered) %*% (centered * weights)
}

#' Construct an importance-sampling approximation to the ZIZ posterior.
#'
#' @param x An input object or numeric vector required by the helper.
#' @param prior A prior object created by `makePrior()`.
#' @param shape1 First beta-prior shape parameter for the zero-inflation probability.
#' @param shape2 Second beta-prior shape parameter for the zero-inflation probability.
#' @param nSamples Number of importance or posterior samples.
#' @param proposalScale Positive multiplier applied to the Laplace covariance for the importance proposal.
#' @param seed Optional random-number seed.
#' @param start Starting values for optimisation or fitting.
#' @param laplace Optional precomputed Laplace approximation used to construct the proposal.
#' @return A list describing weighted samples, posterior moments, proposal, and diagnostics.
#' @keywords internal
#' @noRd
makeZizPosteriorImportance = function(x,
                                      prior,
                                      shape1 = 1,
                                      shape2 = 1,
                                      nSamples = 5000,
                                      proposalScale = 2,
                                      seed = NULL,
                                      start = c(pi = 0.5, shape = 2),
                                      laplace = NULL) {
  validateBayesPrior(prior)

  if (!is.numeric(shape1) || length(shape1) != 1L || !is.finite(shape1) || shape1 <= 0) {
    stop("shape1 must be a positive finite number")
  }

  if (!is.numeric(shape2) || length(shape2) != 1L || !is.finite(shape2) || shape2 <= 0) {
    stop("shape2 must be a positive finite number")
  }

  nSamples = as.integer(nSamples)
  if (!is.finite(nSamples) || nSamples < 100L) {
    stop("nSamples must be at least 100")
  }

  if (!is.numeric(proposalScale) || length(proposalScale) != 1L ||
      !is.finite(proposalScale) || proposalScale <= 0) {
    stop("proposalScale must be a positive finite number")
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("seed must be NULL or a finite numeric value")
    }

    set.seed(as.integer(seed))
  }

  if (is.null(laplace)) {
    laplace = makeZizPosteriorLaplace(
      x = x,
      prior = prior,
      shape1 = shape1,
      shape2 = shape2,
      start = start
    )
  }

  if (!is.list(laplace) || is.null(laplace$modeWorking) || is.null(laplace$covarianceWorking)) {
    stop("laplace must be a zero-inflated zeta Laplace approximation")
  }

  proposalMean = unname(laplace$modeWorking)
  names(proposalMean) = c("eta", "tau")
  # Inflate the local Laplace covariance so the proposal has heavier practical
  # coverage of posterior tails; importance weights correct for the mismatch.
  proposalCovariance = laplace$covarianceWorking * proposalScale
  dimnames(proposalCovariance) = list(c("eta", "tau"), c("eta", "tau"))

  workingSamples = makeZizProposalDraws(
    mean = proposalMean,
    covariance = proposalCovariance,
    n = nSamples
  )

  obsData = zizObservationData(x)
  counts = x$data$rn
  logPosterior = apply(
    workingSamples,
    1L,
    zizWorkingLogPosterior,
    obsData = obsData,
    counts = counts,
    prior = prior,
    shape1 = shape1,
    shape2 = shape2
  )
  logProposal = zizProposalLogDensity(
    working = workingSamples,
    mean = proposalMean,
    covariance = proposalCovariance
  )

  # Importance weights are posterior/proposal ratios. Work on the log scale
  # first, then subtract the maximum finite log weight before exponentiating
  # to avoid numerical overflow and underflow.
  logWeights = logPosterior - logProposal
  finiteWeights = is.finite(logWeights)
  if (!any(finiteWeights)) {
    stop("importance sampling produced no finite posterior weights")
  }

  logWeightScale = max(logWeights[finiteWeights])
  scaledWeights = exp(logWeights - logWeightScale)
  scaledWeights[!is.finite(scaledWeights)] = 0
  weightSum = sum(scaledWeights)

  if (!is.finite(weightSum) || weightSum <= 0) {
    stop("importance sampling weights could not be normalized")
  }

  # Normalised weights sum to one and can therefore be used directly for
  # posterior moments, probability summaries, and effective sample size.
  weights = scaledWeights / weightSum
  thetaSamples = t(apply(workingSamples, 1L, zizWorkingToTheta))
  colnames(thetaSamples) = c("pi", "shape")
  posteriorMean = colSums(thetaSamples * weights)
  names(posteriorMean) = c("pi", "shape")
  posteriorCovariance = weightedCovariance(thetaSamples, weights, posteriorMean)
  dimnames(posteriorCovariance) = list(c("pi", "shape"), c("pi", "shape"))

  effectiveSampleSize = 1 / sum(weights^2)
  maxWeight = max(weights)

  weightedSamples = data.frame(
    pi = thetaSamples[, "pi"],
    shape = thetaSamples[, "shape"],
    eta = workingSamples[, "eta"],
    tau = workingSamples[, "tau"],
    logPosterior = logPosterior,
    logProposal = logProposal,
    weight = weights
  )

  list(
    samples = weightedSamples,
    mean = posteriorMean,
    varCov = posteriorCovariance,
    laplace = laplace,
    proposalMean = proposalMean,
    proposalCovariance = proposalCovariance,
    diagnostics = list(
      effectiveSampleSize = unname(effectiveSampleSize),
      maxWeight = unname(maxWeight),
      nSamples = nSamples,
      proposalScale = proposalScale
    ),
    posteriorMethod = "importance"
  )
}
