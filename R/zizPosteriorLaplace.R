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

  logValue = zizLogLikelihood(obsData, counts, pi, shape) +
    dbeta(pi, shape1, shape2, log = TRUE) +
    prior$logd(shape) +
    zizWorkingLogJacobian(working)

  if (!is.finite(logValue)) {
    return(-Inf)
  }

  unname(logValue)
}

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

fitZIDistBayesLaplace = function(x,
                                 nterms = 10,
                                 prior = makePrior(),
                                 shape1 = 1,
                                 shape2 = 1,
                                 start = c(pi = 0.5, shape = 2),
                                 ...) {
  nvals = 1:nterms
  approximation = makeZizPosteriorLaplace(
    x = x,
    prior = prior,
    shape1 = shape1,
    shape2 = shape2,
    start = start
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
    model = "ziz",
    method = "bayes",
    posteriorMethod = "laplace"
  )

  class(result) = "psFit"
  result
}
