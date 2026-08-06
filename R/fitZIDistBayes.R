#' Fit a zero-inflated zeta model by Metropolis-Hastings sampling
#'
#' @param x An object of class `psData`.
#' @param nterms Number of fitted probabilities to retain.
#' @param prior A shape prior returned by [makePrior()].
#' @param theta0 Initial values for `pi` and `shape`.
#' @param shape1 First shape parameter of the beta prior for `pi`.
#' @param shape2 Second shape parameter of the beta prior for `pi`.
#' @param nIter Number of retained MCMC iterations.
#' @param nBurnIn Number of burn-in iterations.
#' @param silent Logical; suppress the progress bar when `TRUE`.
#' @param seed Optional integer seed.
#' @param ... Retained for backward compatibility.
#'
#' @return An object of class `psFit`.
#'
#' @keywords internal
#' @importFrom methods is
#' @importFrom stats cov dbeta rbeta runif
#' @importFrom utils setTxtProgressBar txtProgressBar
fitZIDistBayes = function(x,
                          nterms = 10,
                          prior = makePrior(),
                          theta0 = c(0.5, 2),
                          shape1 = 1,
                          shape2 = 1,
                          nIter = 1e4,
                          nBurnIn = 1e3,
                          silent = TRUE,
                          seed = NULL,
                          ...) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData")
  }

  dotargs = list(...)
  if (missing(prior) && any(c("a", "b") %in% names(dotargs))) {
    a = if ("a" %in% names(dotargs)) dotargs$a else -2
    b = if ("b" %in% names(dotargs)) dotargs$b else 2

    if (!is.numeric(a) || length(a) != 1L || !is.finite(a) ||
        !is.numeric(b) || length(b) != 1L || !is.finite(b) || b <= a) {
      stop("Legacy MCMC bounds must be finite numbers with b greater than a")
    }

    prior = makePrior(
      family = "loguniform",
      range = 1 + exp(c(a, b))
    )
  }

  validateBayesPrior(prior)

  nIter = as.integer(nIter)
  nBurnIn = as.integer(nBurnIn)

  if (!is.finite(nIter) || nIter <= 0L) {
    stop("nIter must be greater than zero")
  }

  if (!is.finite(nBurnIn) || nBurnIn <= 0L) {
    stop("nBurnIn must be greater than zero")
  }

  if (nIter < 1000L) {
    warning(
      "The number of retained MCMC samples should usually be 1000 or higher.",
      call. = FALSE
    )
  }

  if (!is.numeric(theta0) || length(theta0) != 2L || any(!is.finite(theta0))) {
    stop("theta0 must be a finite numeric vector of length two")
  }

  if (theta0[1] <= 0 || theta0[1] >= 1) {
    stop("The starting value for pi must be in (0, 1)")
  }

  validateZetaShape(theta0[2], "theta0 shape")

  if (!inRange(theta0[2], prior$range)) {
    stop("The starting shape value must lie strictly inside the prior range")
  }

  if (!is.numeric(shape1) || length(shape1) != 1L ||
      !is.finite(shape1) || shape1 <= 0) {
    stop("shape1 must be a positive finite number")
  }

  if (!is.numeric(shape2) || length(shape2) != 1L ||
      !is.finite(shape2) || shape2 <= 0) {
    stop("shape2 must be a positive finite number")
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("seed must be NULL or one finite number")
    }
    set.seed(as.integer(seed))
  }

  obsData = zizObservationData(x)
  counts = x$data$rn
  nTotal = nIter + nBurnIn
  currentPi = unname(theta0[1])
  currentShape = unname(theta0[2])
  chain = matrix(
    NA_real_,
    nrow = nIter,
    ncol = 2L,
    dimnames = list(NULL, c("pi", "shape"))
  )
  updateShape = sample(c(TRUE, FALSE), nTotal, replace = TRUE)
  logUniforms = log(runif(nTotal))
  acceptedShape = 0L
  acceptedPi = 0L
  proposedShape = 0L
  proposedPi = 0L

  if (!silent) {
    progressBar = txtProgressBar(
      min = 0,
      max = nTotal,
      initial = 0,
      style = 3
    )
    on.exit(close(progressBar), add = TRUE)
  }

  for (iteration in seq_len(nTotal)) {
    if (updateShape[iteration]) {
      proposedShape = proposedShape + 1L
      candidateShape = runif(1L, prior$range[1], prior$range[2])

      logAcceptance =
        zizLogLikelihood(obsData, counts, currentPi, candidateShape) -
        zizLogLikelihood(obsData, counts, currentPi, currentShape) +
        prior$logd(candidateShape) - prior$logd(currentShape)

      if (is.finite(logAcceptance) &&
          (logAcceptance >= 0 || logUniforms[iteration] < logAcceptance)) {
        currentShape = candidateShape
        acceptedShape = acceptedShape + 1L
      }
    } else {
      proposedPi = proposedPi + 1L
      candidatePi = rbeta(1L, shape1 = shape1, shape2 = shape2)

      logTargetRatio =
        zizLogLikelihood(obsData, counts, candidatePi, currentShape) -
        zizLogLikelihood(obsData, counts, currentPi, currentShape) +
        dbeta(candidatePi, shape1, shape2, log = TRUE) -
        dbeta(currentPi, shape1, shape2, log = TRUE)
      logProposalRatio =
        dbeta(currentPi, shape1, shape2, log = TRUE) -
        dbeta(candidatePi, shape1, shape2, log = TRUE)
      logAcceptance = logTargetRatio + logProposalRatio

      if (is.finite(logAcceptance) &&
          (logAcceptance >= 0 || logUniforms[iteration] < logAcceptance)) {
        currentPi = candidatePi
        acceptedPi = acceptedPi + 1L
      }
    }

    if (iteration > nBurnIn) {
      chain[iteration - nBurnIn, ] = c(currentPi, currentShape)
    }

    if (!silent) {
      setTxtProgressBar(progressBar, iteration)
    }
  }

  chain = as.data.frame(chain)
  par = c(pi = mean(chain$pi), shape = mean(chain$shape))
  nvals = seq_len(nterms)
  fitted = (1 - par[["pi"]]) * dzetaStandard(nvals, shape = par[["shape"]])
  fitted[nvals == 1L] = fitted[nvals == 1L] + par[["pi"]]
  names(fitted) = if (x$type == "P") {
    paste0("P", nvals - 1L)
  } else {
    paste0("S", nvals)
  }

  result = list(
    psData = x,
    fit = list(
      par = par,
      acceptance = c(
        pi = if (proposedPi > 0L) acceptedPi / proposedPi else NA_real_,
        shape = if (proposedShape > 0L) acceptedShape / proposedShape else NA_real_
      )
    ),
    pi = unname(par[["pi"]]),
    shape = unname(par[["shape"]]),
    var.cov = cov(chain),
    fitted = fitted,
    chain = chain,
    model = "ziz",
    method = "bayes",
    posteriorMethod = "mcmc"
  )

  class(result) = "psFit"
  result
}
