#' Fit an MCMC posterior for the zero-inflated zeta model.
#'
#' The sampler preserves the established ZIZ Metropolis-Hastings calculation
#' and RNG ordering while returning the common MCMC posterior representation.
#'
#' @param model A `zizModel` descriptor.
#' @param engine An `mcmcPosteriorEngine` object.
#' @param x An object of class `psData`.
#' @param prior A prior object created by `makePrior()`.
#' @param theta0 Initial values for `pi` and `shape`.
#' @param shape1 First shape parameter of the beta prior for `pi`.
#' @param shape2 Second shape parameter of the beta prior for `pi`.
#' @param nIter Number of retained MCMC iterations.
#' @param nBurnIn Number of burn-in iterations.
#' @param silent Logical; suppress the progress bar when `TRUE`.
#' @param seed Optional integer seed.
#' @param ... Additional MCMC controls reserved for future use.
#' @return An `mcmcPosteriorRepresentation` object.
#' @keywords internal
#' @exportS3Method fitMcmcPosteriorModel zizModel
#' @noRd
#' @importFrom stats cov dbeta rbeta runif
#' @importFrom utils setTxtProgressBar txtProgressBar
fitMcmcPosteriorModel.zizModel = function(model,
                                          engine,
                                          x,
                                          prior,
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
  if (!inherits(prior, "psPrior")) {
    stop("prior must be an object of class psPrior")
  }

  validateZetaPriorRange(prior$range)

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

  obsData = modelObservationData(model, x)
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

  # Preserve the legacy RNG sequence by drawing the update indicators and
  # acceptance uniforms before drawing iteration-specific proposals.
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
        modelLogLikelihood(
          model,
          parameters = list(pi = currentPi, shape = candidateShape),
          data = x
        ) -
        modelLogLikelihood(
          model,
          parameters = list(pi = currentPi, shape = currentShape),
          data = x
        ) +
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
        modelLogLikelihood(
          model,
          parameters = list(pi = candidatePi, shape = currentShape),
          data = x
        ) -
        modelLogLikelihood(
          model,
          parameters = list(pi = currentPi, shape = currentShape),
          data = x
        ) +
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
  posteriorMean = c(pi = mean(chain$pi), shape = mean(chain$shape))
  posteriorVariance = cov(chain)
  acceptance = c(
    pi = if (proposedPi > 0L) acceptedPi / proposedPi else NA_real_,
    shape = if (proposedShape > 0L) acceptedShape / proposedShape else NA_real_
  )

  newPsPosteriorRepresentation(
    engine = engine,
    value = list(
      chain = chain,
      mean = posteriorMean,
      variance = posteriorVariance
    ),
    metadata = list(
      model = model$model,
      nIter = nIter,
      nBurnIn = nBurnIn,
      acceptance = acceptance,
      shape1 = shape1,
      shape2 = shape2,
      seed = seed
    )
  )
}
