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
#' @param level Equal-tailed credible interval level for posterior probabilities.
#' @param ... Retained for backward compatibility.
#'
#' @return An object of class `psFit`.
#'
#' @keywords internal
#' @importFrom methods is
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
                          level = 0.95,
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

  model = zizModel()
  engine = mcmcPosteriorEngine()
  representation = fitPosterior(
    engine = engine,
    model = model,
    x = x,
    prior = prior,
    theta0 = theta0,
    shape1 = shape1,
    shape2 = shape2,
    nIter = nIter,
    nBurnIn = nBurnIn,
    silent = silent,
    seed = seed
  )

  diagnostics = posteriorDiagnostics(engine, representation)
  chain = representation$value$chain

  finaliseBayesianPsFit(
    model = model,
    engine = engine,
    representation = representation,
    x = x,
    nterms = nterms,
    level = level,
    fit = list(
      par = posteriorPointEstimate(engine, model, representation),
      acceptance = diagnostics$acceptance
    ),
    legacyFields = list(
      var.cov = representation$value$variance,
      chain = chain
    ),
    posteriorDiagnosticsValue = diagnostics$acceptance,
    useEngineDiagnostics = FALSE
  )
}
