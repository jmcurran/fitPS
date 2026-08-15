#' Fit a logarithmic distribution to forensic data
#'
#' Fits the logarithmic-series distribution to P- or S-survey data using
#' maximum likelihood or the common fitPS Bayesian posterior architecture.
#' The logarithmic model has probability mass function
#' \deqn{p(k) = -\frac{\pi^k}{k\log(1-\pi)}, \quad 0 < \pi < 1.}
#'
#' @param x An object of class `psData`.
#' @param nterms Number of fitted P/S probability terms to retain.
#' @param method Fitting method, either `"mle"` or `"bayes"`.
#' @param prior Optional `psPrior` used for Bayesian fitting. When omitted, a
#'   uniform prior on `(0.001, 0.999)` is used.
#' @param bayesOptions Optional Bayesian controls. The logarithmic model supports
#'   `"numerical"` and `"mcmc"` posterior engines.
#' @param ... Additional fitting controls. For MLE, `start` is an optional
#'   initial value in `(0, 1)`. For MCMC, `pi0`, `nIter`, `nBurnIn`, and
#'   `silent` are passed to the logarithmic MCMC model method.
#' @return An object of class `psFit`.
#' @importFrom stats optim
#' @importFrom VGAM dlog
#' @export
#'
#' @examples
#' data(Psurveys)
#' fit = fitlogDist(Psurveys$roux)
#' fit
fitlogDist = function(x,
                       nterms = 10,
                       method = c("mle", "bayes"),
                       prior,
                       bayesOptions = NULL,
                       ...) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData")
  }
  modelObservationData(logarithmicModel(), x)

  method = match.arg(method)
  model = logarithmicModel()
  nvals = posteriorProbabilityIndices(x$type, nterms)

  if (identical(method, "mle")) {
    dotargs = list(...)
    start = if ("start" %in% names(dotargs)) {
      dotargs$start
    } else if ("pi" %in% names(dotargs)) {
      dotargs$pi
    } else {
      0.5
    }
    validateLogarithmicPi(start, "start")

    objective = function(pi) {
      -modelLogLikelihood(model, parameters = list(pi = pi), data = x)
    }
    fit = optim(
      par = start,
      fn = objective,
      method = "L-BFGS-B",
      lower = sqrt(.Machine$double.eps),
      upper = 1 - sqrt(.Machine$double.eps),
      hessian = TRUE
    )
    pi = unname(fit$par[1L])
    fittedMatrix = modelProbabilities(
      model,
      parameters = list(pi = pi),
      n = nvals,
      type = x$type
    )
    fitted = as.numeric(fittedMatrix[1L, ])
    names(fitted) = colnames(fittedMatrix)

    result = list(
      psData = x,
      fit = fit,
      pi = pi,
      var.pi = 1 / fit$hessian[1L, 1L],
      fitted = fitted,
      model = model$model,
      modelObject = model,
      method = "mle"
    )
    class(result) = "psFit"
    return(result)
  }

  options = if (missing(prior)) {
    if (is.null(bayesOptions)) {
      bayesOptions = list()
    }
    if (is.null(bayesOptions$prior)) {
      bayesOptions$prior = makePrior(
        family = "uniform",
        range = c(0.001, 0.999)
      )
    }
    normaliseBayesOptions(bayesOptions = bayesOptions)
  } else {
    normaliseBayesOptions(bayesOptions = bayesOptions, prior = prior)
  }

  engine = posteriorEngine(options$posteriorMethod)
  validateEngineModelPair(engine, model)
  validateLogarithmicPriorRange(options$prior$range)

  result = fitBayesianModel(
    model = model,
    posteriorMethod = options$posteriorMethod,
    x = x,
    prior = options$prior,
    nterms = nterms,
    ...
  )
  result$bayesOptions = options
  result
}

#' @rdname fitlogDist
#' @export
fitLogdist = fitlogDist

#' @rdname fitlogDist
#' @export
fitlogdist = fitlogDist
