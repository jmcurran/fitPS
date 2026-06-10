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

weightedCovariance = function(values, weights, mean) {
  centered = sweep(values, 2L, mean, `-`)
  t(centered) %*% (centered * weights)
}

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

fitZIDistBayesImportance = function(x,
                                    nterms = 10,
                                    prior = makePrior(),
                                    shape1 = 1,
                                    shape2 = 1,
                                    nSamples = 5000,
                                    proposalScale = 2,
                                    seed = NULL,
                                    start = c(pi = 0.5, shape = 2),
                                    ...) {
  nvals = 1:nterms
  approximation = makeZizPosteriorImportance(
    x = x,
    prior = prior,
    shape1 = shape1,
    shape2 = shape2,
    nSamples = nSamples,
    proposalScale = proposalScale,
    seed = seed,
    start = start
  )

  par = approximation$mean
  fitted = (1 - par[["pi"]]) * dzetaStandard(nvals, shape = par[["shape"]])
  fitted[nvals == 1] = fitted[nvals == 1] + par[["pi"]]
  names(fitted) = if (x$type == "P") {
    paste0("P", nvals - 1)
  } else {
    paste0("S", nvals)
  }

  result = list(
    psData = x,
    fit = list(par = par),
    pi = unname(par[["pi"]]),
    shape = unname(par[["shape"]]),
    var.cov = approximation$varCov,
    fitted = fitted,
    weightedSamples = approximation$samples,
    importance = approximation,
    posteriorDiagnostics = approximation$diagnostics,
    model = "ziz",
    method = "bayes",
    posteriorMethod = "importance"
  )

  class(result) = "psFit"
  result
}
