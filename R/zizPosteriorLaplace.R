#' Evaluate the ZIZ log posterior in unconstrained working coordinates.
#'
#' @param working Numeric vector or matrix of working parameters `(eta, tau)`.
#' @param obsData Positive-integer observation support used by the ZIZ likelihood.
#' @param counts Observation frequencies corresponding to `obsData`.
#' @param prior A prior object created by `makePrior()`.
#' @param shape1 First beta-prior shape parameter for the zero-inflation probability.
#' @param shape2 Second beta-prior shape parameter for the zero-inflation probability.
#' @return A numeric log posterior value on the working scale, or `-Inf` outside the valid region.
#' @keywords internal
#' @noRd
zizWorkingLogPosterior = function(working,
                                  obsData,
                                  counts,
                                  prior,
                                  shape1 = 1,
                                  shape2 = 1) {
  working = unname(working)

  if (!is.numeric(working) || length(working) != 2L || any(!is.finite(working))) {
    return(-Inf)
  }

  theta = zizWorkingToTheta(working)
  pi = unname(theta[["pi"]])
  shape = unname(theta[["shape"]])

  # The target density is expressed in working coordinates, so the natural-
  # scale posterior must include the transformation Jacobian. Omitting this
  # term would approximate a different posterior after reparameterisation.
  logValue = zizLogLikelihood(obsData, counts, pi, shape) +
    dbeta(pi, shape1, shape2, log = TRUE) +
    prior$logd(shape) +
    zizWorkingLogJacobian(working)

  if (!is.finite(logValue)) {
    return(-Inf)
  }

  unname(logValue)
}

#' Construct a Laplace approximation to the ZIZ posterior.
#'
#' @param x An input object or numeric vector required by the helper.
#' @param prior A prior object created by `makePrior()`.
#' @param shape1 First beta-prior shape parameter for the zero-inflation probability.
#' @param shape2 Second beta-prior shape parameter for the zero-inflation probability.
#' @param start Starting values for optimisation or fitting.
#' @return A list describing the posterior mode, Hessian, and covariance approximations.
#' @keywords internal
#' @noRd
makeZizPosteriorLaplace = function(x,
                                   prior,
                                   shape1 = 1,
                                   shape2 = 1,
                                   start = c(pi = 0.5, shape = 2)) {
  validateBayesPrior(prior)

  if (!is.numeric(shape1) || length(shape1) != 1L || !is.finite(shape1) || shape1 <= 0) {
    stop("shape1 must be a positive finite number")
  }

  if (!is.numeric(shape2) || length(shape2) != 1L || !is.finite(shape2) || shape2 <= 0) {
    stop("shape2 must be a positive finite number")
  }

  if (!is.numeric(start) || length(start) != 2L || any(!is.finite(start))) {
    stop("start must be a finite numeric vector of length two")
  }

  if (start[[1L]] <= 0 || start[[1L]] >= 1) {
    stop("start pi must be strictly between 0 and 1")
  }

  validateZetaShape(start[[2L]], "start shape")

  if (!inRange(start[[2L]], prior$range)) {
    start[[2L]] = mean(prior$range)
  }

  obsData = zizObservationData(x)
  counts = x$data$rn
  startWorking = zizThetaToWorking(start)

  objective = function(working) {
    logValue = zizWorkingLogPosterior(
      working = working,
      obsData = obsData,
      counts = counts,
      prior = prior,
      shape1 = shape1,
      shape2 = shape2
    )

    if (!is.finite(logValue)) {
      return(.Machine$double.xmax^0.25)
    }

    -logValue
  }

  workingLower = c(eta = -36, tau = log(prior$range[[1L]] - 1))
  workingUpper = c(eta = 36, tau = log(prior$range[[2L]] - 1))

  fit = optim(
    par = startWorking,
    fn = objective,
    method = "L-BFGS-B",
    lower = workingLower,
    upper = workingUpper,
    hessian = TRUE
  )

  if (!isTRUE(fit$convergence == 0L)) {
    stop("Laplace posterior optimisation did not converge")
  }

  hessian = fit$hessian
  if (!is.matrix(hessian) || any(!is.finite(hessian))) {
    stop("Laplace posterior Hessian is not finite")
  }

  # At the posterior mode, the inverse Hessian of the negative log posterior
  # gives the local Gaussian covariance used by the Laplace approximation.
  covarianceWorking = tryCatch(
    solve(hessian),
    error = function(e) {
      NULL
    }
  )

  if (is.null(covarianceWorking) || any(!is.finite(covarianceWorking))) {
    stop("Laplace posterior covariance could not be computed")
  }

  eigenValues = eigen(covarianceWorking, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigenValues <= 0)) {
    stop("Laplace posterior covariance is not positive definite")
  }

  modeWorking = unname(fit$par)
  names(modeWorking) = c("eta", "tau")
  mode = zizWorkingToTheta(modeWorking)

  jacobianTheta = matrix(
    c(
      mode[["pi"]] * (1 - mode[["pi"]]), 0,
      0, mode[["shape"]] - 1
    ),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("pi", "shape"), c("eta", "tau"))
  )
  # Transform the local working-scale covariance back to (pi, shape) using
  # the first-order delta method.
  covarianceTheta = jacobianTheta %*% covarianceWorking %*% t(jacobianTheta)

  dimnames(covarianceWorking) = list(c("eta", "tau"), c("eta", "tau"))
  dimnames(covarianceTheta) = list(c("pi", "shape"), c("pi", "shape"))

  list(
    mode = mode,
    modeWorking = modeWorking,
    covarianceWorking = covarianceWorking,
    varCov = covarianceTheta,
    hessian = hessian,
    logPosteriorMode = -fit$value,
    convergence = fit$convergence,
    counts = counts,
    obsData = obsData,
    posteriorMethod = "laplace"
  )
}

#' Fit the zero-inflated zeta model using a Laplace posterior approximation.
#'
#' @param x An input object or numeric vector required by the helper.
#' @param nterms Number of fitted P/S probability terms to retain.
#' @param prior A prior object created by `makePrior()`.
#' @param shape1 First beta-prior shape parameter for the zero-inflation probability.
#' @param shape2 Second beta-prior shape parameter for the zero-inflation probability.
#' @param start Starting values for optimisation or fitting.
#' @param nPosteriorDraws Number of draws generated from the Laplace approximation for posterior summaries.
#' @param seed Optional random-number seed.
#' @param level Probability level for intervals or summaries.
#' @param ... Additional arguments passed to the underlying fitting or helper routine.
#' @return A Bayesian ZIZ `psFit` object.
#' @keywords internal
#' @noRd
fitZIDistBayesLaplace = function(x,
                                 nterms = 10,
                                 prior = makePrior(),
                                 shape1 = 1,
                                 shape2 = 1,
                                 start = c(pi = 0.5, shape = 2),
                                 nPosteriorDraws = 5000,
                                 seed = NULL,
                                 level = 0.95,
                                 ...) {
  nvals = 1:nterms
  approximation = makeZizPosteriorLaplace(
    x = x,
    prior = prior,
    shape1 = shape1,
    shape2 = shape2,
    start = start
  )

  nPosteriorDraws = as.integer(nPosteriorDraws)
  if (!is.finite(nPosteriorDraws) || nPosteriorDraws < 100L) {
    stop("nPosteriorDraws must be at least 100")
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("seed must be NULL or a finite numeric value")
    }
    set.seed(as.integer(seed))
  }

  workingDraws = makeZizProposalDraws(
    mean = approximation$modeWorking,
    covariance = approximation$covarianceWorking,
    n = nPosteriorDraws
  )
  thetaDraws = t(apply(workingDraws, 1L, zizWorkingToTheta))
  posteriorProbs = summariseZizSampleProbabilities(
    pi = thetaDraws[, "pi"],
    shape = thetaDraws[, "shape"],
    type = x$type,
    nterms = nterms,
    level = level,
    posteriorMethod = "laplace"
  )

  par = approximation$mode
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
    laplace = approximation,
    posteriorProbs = posteriorProbs,
    model = "ziz",
    method = "bayes",
    posteriorMethod = "laplace"
  )

  result = attachZizPosterior(
    result = result,
    probabilities = posteriorProbs,
    representation = approximation,
    diagnostics = NULL,
    level = level
  )

  class(result) = "psFit"
  result
}
