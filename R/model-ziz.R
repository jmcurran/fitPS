# Zero/one-inflated zeta model implementation
#
# Distribution-specific ZIZ fitting, probability, posterior, simulation, and
# model-dispatch code is intentionally consolidated in this file.

# ---- fitZIDist.R ----
#' Fit a Zero-Inflated Zeta Distribution to Forensic Data
#'
#' This function uses maximum likelihood estimation (MLE) or Bayesian estimation
#' to estimate the mixing parameter and the shape parameter of a
#' zero-inflated zeta distribution from a set of observed counts for either the
#' number of groups/sources of forensically interesting material (mostly glass
#' or paint) recovered from clothing, or the number of fragments/particles in
#' each group. This, in turn, allows the estimation of the P and S
#' probabilities, as described by Evett and Buckleton (1990), which are used in
#' computing the likelihood ratio (LR) for activity level propositions. The data
#' arise from clothing surveys. The zero-inflated zeta distribution has
#' probability mass function
#' \deqn{p(k) = \begin{cases}
#' \pi + \frac{(1-\pi)}{\zeta(s)}&,k=0, \\
#' \frac{(1-\pi)k^{-s}}{\zeta(s)}&,k=1,2,\ldots
#' \end{cases}
#' }{
#' p(0) = pi + (1-pi)/zeta(s),k = 0 and p(k) = (1-pi)k^-s/zeta(s),k>1
#' }
#' where \eqn{\zeta(s)}{zeta(s)} is the Riemann zeta function.
#'
#' @details The function returns an object of class \code{psFit}. Core fitted
#'   values are stored at the top level, while Bayesian uncertainty is stored
#'   in the attached \code{psPosterior} object. Important components include:
#' \describe{
#' \item{\code{psData}}{ -- an object of class \code{psData}--see \code{\link{readData}},}
#' \item{\code{fit}}{ -- method-specific point-estimate information,}
#' \item{\code{pi}}{ - the maximum likelihood estimate, or the posterior mean, of the mixing parameter,}
#' \item{\code{shape}}{ -- the maximum likelihood estimate, or the posterior mean, of the shape parameter,}
#' \item{\code{var.cov}}{ -- for maximum-likelihood fits, the estimated variance-covariance matrix for the parameters,}
#' \item{\code{fitted}}{ -- a named \code{vector} containing the first \code{nterms} of the fitted distribution.}
#' \item{\code{model}}{ -- set to \code{"ziz"} for this model,}
#' \item{\code{method}}{ -- the method of estimation used, either \code{"mle"} or \code{"bayes"},}
#' \item{\code{posterior}}{ -- for Bayesian fits, an object of class \code{psPosterior} containing posterior parameter summaries, posterior probability summaries, and the engine-specific posterior representation,}
#' }
#'
#' The output can be used in a variety of ways. If the interest is just in the
#' mixing and shape parameter estimates, then the \code{pi} and \code{shape}
#' members of the \code{psFit} object contain this information. It is also
#' displayed along with a number of fitted probabilities by the
#' \code{\link{print.psFit}} method. The fitted object can also be plotted using
#' the plot method \code{\link{plot.psFit}}, and to create a probability
#' function with \code{\link{probfun}}. The \code{shape} value stored in
#' the fitted object is the zeta distribution shape parameter and must satisfy
#' \code{shape > 1}.
#'
#' This function implements both maximum likelihood estimation (MLE) and
#' Bayesian estimation. Both modes of estimation require additional information
#' such as starting values and parameters for priors. Please read the
#' documentation for the \code{...} argument closely because it explains what
#' you can change and what the default values are.
#'
#' Bayesian zero-inflated zeta estimation is selected with
#' \code{method = "bayes"}. The posterior approximation is selected with
#' \code{bayesOptions$posteriorMethod}. The default, \code{"numerical"}, uses
#' deterministic two-dimensional grid integration over \code{pi} and
#' \code{shape}. The legacy Metropolis-Hastings sampler remains available with
#' \code{bayesOptions = list(posteriorMethod = "mcmc")}. The prior for the
#' mixing proportion is Beta(shape1, shape2), and the shape prior is supplied
#' by \code{bayesOptions$prior} or \code{prior}. If no shape prior is supplied,
#' \code{makePrior()} is used. Bayesian fits store a \code{psPosterior}
#' object in \code{fit$posterior}; use \code{posteriorProbs()} for posterior
#' means and credible intervals of the derived P or S probabilities. The four
#' posterior engines are \code{"numerical"}, \code{"mcmc"},
#' \code{"laplace"}, and \code{"importance"}.

#'
#' @seealso \code{\link{plot.psFit}}, \code{\link{print.psFit}},
#'   \code{\link{probfun}}, \code{\link{posteriorProbs}},
#'   \code{\link{bootstrapProbs}}.
#'
#' @references Evett, I. W. and Buckleton, J. S., "The interpretation of glass
#'   evidence. A practical approach", Journal of the Forensic Science Society
#'   1990: 30(4): 215--223.
#'
#' @param x an object of type \code{psData}, usually obtained from
#'   \code{\link{readData}}.
#' @param nterms the number of terms to compute the probability distribution
#'   for.
#' @param method primary fitting method. Use \code{"mle"} for maximum
#'   likelihood estimation or \code{"bayes"} for Bayesian estimation. Legacy
#'   Bayesian aliases \code{"integrate"}, \code{"numerical"},
#'   \code{"mcmc"}, \code{"laplace"}, and \code{"importance"} are accepted
#'   with a deprecation warning and translated to \code{method = "bayes"}
#'   with the corresponding \code{bayesOptions$posteriorMethod}.
#' @param prior optional prior object used by Bayesian posterior approximation
#'   methods where applicable. This is retained for consistency with
#'   \code{fitDist()}; new code should usually pass priors through
#'   \code{bayesOptions}.
#' @param bayesOptions optional list controlling Bayesian fitting. The
#'   \code{posteriorMethod} element selects \code{"numerical"},
#'   \code{"mcmc"}, \code{"laplace"}, or \code{"importance"}. The
#'   default is \code{"numerical"}. The \code{prior} element may contain
#'   a prior object returned by \code{makePrior()}.
#' @param ... other arguments that control the estimation methods. If
#'   \code{method == "mle"}, then the user can provide an optional argument
#'   \code{start} which is the starting value for the numerical optimisation. If
#'   this is not provided, then \code{start = c(0.5, 2)} by default. If you specify your
#'   own starting value, keep the mixing parameter greater than 0.5 and use
#'   \code{shape > 1}.
#'
#'   If \code{method == "bayes"}, engine-specific controls can be supplied
#'   through \code{...}. Common MCMC controls include:
#' \describe{
#'   \item{\code{theta0}}{ -- The initial values of the mixing parameter and shape parameter. The default is \code{c(0.5, 2)}. }
#'   \item{\code{shape1}}{ -- The first shape parameter for the beta prior on the mixing distribution, Beta(shape1, shape2). The default is 1. }
#'   \item{\code{shape2}}{ -- The second shape parameter for the beta prior on the mixing distribution, Beta(shape1, shape2). The default is 1. }
#'   \item{\code{nIter}}{ -- The number of samples to save from the chain. Must be greater than zero, and ideally greater than 1000. }
#'   \item{\code{nBurnIn}}{ -- The number of samples to discard from the chain. Must be greater than zero. **NOTE**: the sampler runs for \code{nIter + nBurnIn} iterations,
#'   so you do not need to factor this number into your number of samples, \code{nIter}. }
#'   \item{\code{silent}}{ -- A logical variable which allows the user to get a progress bar if they want. \code{TRUE} by default. }
#' }
#'
#' @importFrom stats optim runif
#' @importFrom VGAM dzeta
#'
#' @section Deprecated:
#' `fitZIDist()` is retained as a compatibility wrapper. New code should use
#' `fit(x, model = zizModel())`.
#'
#' @return an object of class \code{psFit}--see Details.
#' @export
#'
#' @examples
#' data(Psurveys)
#' roux = Psurveys$roux
#' fit = fitZIDist(roux)
#' fit
fitZIDist = function(x,
                      nterms = 10,
                      method = c(
                        "mle", "bayes", "integrate", "numerical", "mcmc",
                        "laplace", "importance"
                      ),
                      prior,
                      bayesOptions = NULL,
                      ...) {
  signalDeprecatedFitter(
    old = "fitZIDist",
    replacement = "fit(x, model = zizModel())"
  )

  methodInfo = normaliseBayesMethod(method, bayesOptions = bayesOptions)
  args = list(
    x = x,
    model = zizModel(),
    nterms = nterms,
    method = methodInfo$method,
    bayesOptions = methodInfo$bayesOptions,
    ...
  )
  if (!missing(prior)) {
    args$prior = prior
  }

  do.call(fit, args)
}

#' Internal zero-inflated zeta fitting implementation.
#'
#' @inheritParams fitZIDist
#' @return An object of class `psFit`.
#' @keywords internal
#' @noRd
fitZIDistImpl = function(x, nterms = 10,
                     method = c("mle", "bayes", "integrate", "numerical", "mcmc", "laplace", "importance"),
                     prior,
                     bayesOptions = NULL,
                     ...){
  nvals = 1:nterms
  if(!is(x, "psData")){
    stop("x must be an object of class psData")
  }

  if(length(x$data$n) < 2){
    if(x$type == "S"){
      stop("There has to be at least one value higher than 1")
    }else{
      stop("There has to be at least one value higher than 0")
    }
  }

  obsData = if(x$type == 'P'){ ## the main difference is that the values need 1 added
    x$data$n + 1
  }else{
    x$data$n
  }

  methodInfo = normaliseBayesMethod(method, bayesOptions = bayesOptions)
  method = methodInfo$method
  bayesOptions = methodInfo$bayesOptions

  model = zizModel()

  if(method == "mle"){

    dotargs = list(...)
    if("start" %in% names(dotargs)){
      start = dotargs$start
    }else if("shape" %in% names(dotargs)){
      start = dotargs$shape
    }else{
      start = c(0.5, 2)
    }

    if(start[1] <= 0 || start[1] >= 1){
      stop("The starting value for pi must be in (0, 1)")
    }

    validateZetaShape(start[2], "start shape")


    y = rep(obsData, x$data$rn)

    d.one.inflated.zeta = function(x, shape, p, log = FALSE){
      rval = (1 - p) * dzetaStandard(x, shape = shape)
      rval[x == 1] = rval[x == 1] + p

      if(log){
        return(log(rval))
      }
      return(rval)
    }

    logLik = function(params){
      p = params[1]
      shape = params[2]

      validateZetaShape(shape)
      rval = (1 - p) * dzetaStandard(obsData, shape = shape)
      rval[obsData == 1] = rval[obsData == 1] + p

      r = -sum(x$data$rn * log(rval))
      if(is.infinite(r) || is.nan(r)){
        stop(sprintf("Infinite log-likelihod: pi = %6.4E shape = %6.4f\n", p, shape))
      }
      return(r)
    }

    # fit = nlminb(start = start,
    #              objective = logLik,
    #              lower = 1)

    fit = optim(par = start,
                fn = logLik,
                method = "L-BFGS-B",
                lower = c(sqrt(.Machine$double.eps), 1 + sqrt(.Machine$double.eps)),
                upper  = c(1 - .Machine$double.eps, Inf),
                hessian = TRUE)

    fitted = d.one.inflated.zeta(nvals, shape = fit$par[2], p = fit$par[1])
    names(fitted) = if(x$type == 'P'){
      paste0("P", nvals - 1)
    }else{
      paste0("S", nvals)
    }

    result = list(
      psData = x,
      fit = fit,
      pi = fit$par[1],
      shape =  fit$par[2],
      var.cov = solve(fit$hessian),
      fitted = fitted,
      model = model$model,
      modelObject = model,
      method = "mle"
    )


    class(result) = "psFit"

    return(result)
  } else { ## method == "bayes"
    options = if (missing(prior)) {
      normaliseBayesOptions(bayesOptions = bayesOptions)
    } else {
      normaliseBayesOptions(bayesOptions = bayesOptions, prior = prior)
    }

    engine = posteriorEngine(options$posteriorMethod)
    validateEngineModelPair(engine, model)

    result = fitBayesianModel(
      model = model,
      posteriorMethod = options$posteriorMethod,
      x = x,
      prior = options$prior,
      nterms = nterms,
      ...
    )
    result$bayesOptions = options
    return(result)
  }
}

#' @rdname fitZIDist
#' @export
fitZIdist = fitZIDist

#' @rdname fitZIDist
#' @export
fitzidist = fitZIDist

#' Evaluate the legacy zero-inflated zeta log likelihood used by profile-likelihood fitting.
#'
#' @param y Observed values used by the legacy likelihood helper.
#' @param theta Numeric vector containing natural ZIZ parameters `(pi, shape)`.
#' @return A numeric log-likelihood value.
#' @keywords internal
#' @noRd
zi.loglik = function(y, theta){
  p = theta[1]
  shape = theta[2]
  rval = (1 - p) * dzetaStandard(y$n, shape = shape)
  rval[y$n == 1] = rval[y$n == 1] + p
  sum(y$rn *log(rval))
}


#' Fit the zero-inflated zeta model by profile likelihood.
#'
#' @param x An input object or numeric vector required by the helper.
#' @param nterms Number of fitted P/S probability terms to retain.
#' @param start Starting values for optimisation or fitting.
#' @param lambda Penalty or tuning value used by the legacy profile-likelihood fit.
#' @param ... Additional arguments passed to the underlying fitting or helper routine.
#' @return A zero-inflated zeta `psFit` object.
#' @keywords internal
#' @noRd
fitZIDistPL = function(x, nterms = 10,
                       start = c(0.5, 2),
                       lambda = 0.1,
                     ...){
  nvals = 1:nterms
  if(!is(x, "psData")){
    stop("x must be an object of class psData")
  }

  if(start[1] <= 0 || start[1] >= 1){
    stop("The starting value for pi must be in (0, 1)")
  }

  validateZetaShape(start[2], "start shape")

  if(length(x$data$n) < 2){
    if(x$type == "S"){
      stop("There has to be at least one value higher than 1")
    }else{
      stop("There has to be at least one value higher than 0")
    }
  }

  obsData = if(x$type == 'P'){ ## the main difference is that the values need 1 added
    x$data$n + 1
  }else{
    x$data$n
  }

  y = rep(obsData, x$data$rn)

  d.one.inflated.zeta = function(x, shape, p, log = FALSE){
    rval = (1 - p) * dzetaStandard(x, shape = shape)
    rval[x == 1] = rval[x == 1] + p

    if(log){
      return(log(rval))
    }
    return(rval)
  }

  logLik = function(params){
    p = params[1]
    shape = params[2]

    validateZetaShape(shape)
    rval = (1 - p) * dzetaStandard(obsData, shape = shape)
    rval[obsData == 1] = rval[obsData == 1] + p

    r = -(sum(x$data$rn * log(rval)) - lambda * (log(p)  + log(1-p)))
    if(is.infinite(r) || is.nan(r)){
      stop(sprintf("Infinite log-likelihod: pi = %6.4E shape = %6.4f\n", p, shape))
    }
    return(r)
  }

  # fit = nlminb(start = start,
  #              objective = logLik,
  #              lower = 1)

  fit = optim(par = start,
              fn = logLik,
              method = "L-BFGS-B",
              lower = c(sqrt(.Machine$double.eps), 1 + sqrt(.Machine$double.eps)),
              upper  = c(1 - .Machine$double.eps, Inf),
              hessian = TRUE)

  fitted = d.one.inflated.zeta(nvals, shape = fit$par[2], p = fit$par[1])
  names(fitted) = if(x$type == 'P'){
    paste0("P", nvals - 1)
  }else{
    paste0("S", nvals)
  }

  result = list(
    psData = x,
    fit = fit,
    pi = fit$par[1],
    shape =  fit$par[2],
    var.cov = solve(fit$hessian),
    fitted = fitted,
    model = "ziz",
    modelObject = zizModel()
  )


  class(result) = "psFit"

  return(result)
}

# ---- profileLikelihoodZIZ.R ----
#' Construct a profile-likelihood surface and interval information for a zero-inflated zeta fit.
#'
#' @param x An input object or numeric vector required by the helper.
#' @param level Probability level for intervals or summaries.
#' @param grid.Pi Grid of zero-inflation probabilities used for profiling.
#' @param grid.Shape Grid of zeta shape values used for profiling.
#' @param silent Logical; suppress progress output when `TRUE`.
#' @return An internal profile-likelihood result.
#' @keywords internal
#' @noRd
profileLikelihoodZIZ  = function(x, level = 0.95,
                 grid.Pi = seq(0 + .Machine$double.eps, 1 - .Machine$double.eps, length = 100),
                 grid.Shape = seq(1 + sqrt(.Machine$double.eps), x$shape + 4 * sqrt(diag(x$var.cov))[2], by = 0.01),
                 silent = FALSE){
  if(x$model != "ziz"){
    stop("This function is for the ZIZ model only")
  }

  if(any(level <= 0 | level >= 1)){
    stop("The elements of level must be in the interval (0, 1)")
  }

  sx = sqrt(diag(x$var.cov))[2]
  l0 = -x$fit$val
  dataf = x$psData$data
  if(x$psData$type == "P"){
    dataf$n = dataf$n + 1
  }

  if(!silent){
    message("Computing contours. This may take a few seconds.")
  }

  r = lapply(grid.Pi, function(p){
    sapply(grid.Shape, function(si){
      -2 * (zi.loglik(dataf, c(p, si)) - l0)
    })
  })

  r = do.call("rbind", r)

  cr = lapply(level, function(l){
    qstar2 = qchisq(l, 2)
    r1 = r - qstar2
    confRegion = contourLines(grid.Pi, grid.Shape, r1, levels = 0)[[1]]
    confRegion = data.frame(pi = confRegion$x, shape = confRegion$y)
  })

  names(cr) = as.character(level)

  return(cr)
}

# ---- rZIzeta.R ----
#' Generate zero inflated zeta random variates
#'
#' @param n the number of observations.
#' @param pi the mixing parameter for the zero-inflated zeta model---must be in
#'   (0, 1).
#' @param shape the shape parameter for the zero-inflated zeta. Must be greater
#'   than 1.
#' @param offset the zeta distribution returns random variates that are greater
#'   than, or equal to one. If the offset is greater than 0, then the
#'   distribution is anchored on (has minimum value of) \code{1 - offset}.
#'
#' @return a vector of random variates from a zero-inflated zeta model
#'
#' @details Technically this function returns values from the one-inflated zeta
#'   distribution. However, if \code{offset} is greater than zero (and typically
#'   we expect it to be 1), then the minimium random variate value is \code{1 -
#'   offset}. We chose the name "zero-inflated zeta" as more people are familiar
#'   with zero-inflated models.
#'
#'
#' @examples
#' data(Psurveys)
#' roux = Psurveys$roux
#' fit.zi = fitZIDist(roux)
#' x = rZIzeta(n = sum(roux$data$rn), pi = fit.zi$pi, shape = fit.zi$shape)
#' table(x)
#' @export
rZIzeta = function(n, pi = 0.5, shape = 2, offset = 0){

  if(length(pi) > 1 || length(shape) > 1){
    stop("This function does not currently support vector valued inputs for Pi or shape.")
  }

  n = round(n)

  if(n <= 0){
    stop("n must be greater than zero.")
  }

  if(pi <= 0 || pi >= 1){
    stop("Pi must be in (0, 1)")
  }

  validateZetaShape(shape)

  p = runif(n)
  x = p
  x[p < pi] = 1
  n = sum(p >= pi)
  x[p >= pi] = rzetaStandard(n, shape)
  return(x - offset)
}

#' @rdname rZIzeta
#' @export
rzizeta = rZIzeta

#' @rdname rZIzeta
#' @export
rzizeta = rZIzeta
rziz = rZIzeta

# ---- zizNumericalPosterior.R ----
#' Fit a numerical posterior for the zero-inflated zeta model.
#'
#' The ZIZ model retains its existing rectangular two-dimensional grid
#' calculation. The numerical engine wraps that representation in the common
#' shared posterior contract rather than duplicating the integration logic.
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
#' @exportS3Method fitNumericalPosteriorModel zizModel
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

  validateZetaPriorRange(prior$range)

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

# ---- zizMcmcPosterior.R ----
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

# ---- zizPosteriorGrid.R ----
#' Convert P/S data to the positive-integer support used by the ZIZ likelihood.
#'
#' @param x An input object or numeric vector required by the helper.
#' @return A numeric vector on the positive-integer ZIZ support.
#' @keywords internal
#' @noRd
zizObservationData = function(x) {
  psObservationData(x)
}

#' Evaluate the zero-inflated zeta log likelihood for aggregated observations.
#'
#' @param obsData Positive-integer observation support used by the ZIZ likelihood.
#' @param counts Observation frequencies corresponding to `obsData`.
#' @param pi Zero-inflation probability on the natural scale.
#' @param shape Zeta shape parameter on the fitPS scale.
#' @return A numeric log-likelihood value, or `-Inf` outside the parameter space.
#' @keywords internal
#' @noRd
zizLogLikelihood = function(obsData, counts, pi, shape) {
  pi = unname(pi)
  shape = unname(shape)

  if (!is.numeric(pi) || length(pi) != 1L || !is.finite(pi) || pi <= 0 || pi >= 1) {
    return(-Inf)
  }

  if (!is.numeric(shape) || length(shape) != 1L || !is.finite(shape) || shape <= 1) {
    return(-Inf)
  }

  probabilities = (1 - pi) * dzetaStandard(obsData, shape = shape)
  probabilities[obsData == 1] = probabilities[obsData == 1] + pi

  if (any(!is.finite(probabilities)) || any(probabilities <= 0)) {
    return(-Inf)
  }

  sum(counts * log(probabilities))
}

#' Construct and normalise a two-dimensional numerical posterior grid for ZIZ parameters.
#'
#' @param x An input object or numeric vector required by the helper.
#' @param prior A prior object created by `makePrior()`.
#' @param shape1 First beta-prior shape parameter for the zero-inflation probability.
#' @param shape2 Second beta-prior shape parameter for the zero-inflation probability.
#' @param nPiGrid Number of grid points for the zero-inflation probability.
#' @param nShapeGrid Number of grid points for the zeta shape parameter.
#' @return A list describing the normalised joint grid, marginals, moments, and grid spacings.
#' @keywords internal
#' @noRd
makeZizPosteriorGrid = function(x,
                                prior,
                                shape1 = 1,
                                shape2 = 1,
                                nPiGrid = 101,
                                nShapeGrid = 101) {
  validateBayesPrior(prior)

  if (!is.numeric(shape1) || length(shape1) != 1L || !is.finite(shape1) || shape1 <= 0) {
    stop("shape1 must be a positive finite number")
  }

  if (!is.numeric(shape2) || length(shape2) != 1L || !is.finite(shape2) || shape2 <= 0) {
    stop("shape2 must be a positive finite number")
  }

  nPiGrid = as.integer(nPiGrid)
  nShapeGrid = as.integer(nShapeGrid)

  if (!is.finite(nPiGrid) || nPiGrid < 9L) {
    stop("nPiGrid must be at least 9")
  }

  if (!is.finite(nShapeGrid) || nShapeGrid < 9L) {
    stop("nShapeGrid must be at least 9")
  }

  obsData = zizObservationData(x)
  counts = x$data$rn
  # Keep the numerical grid just inside the open support of pi. The clipping
  # avoids log-density failures at exactly 0 or 1 without changing the model.
  piEps = sqrt(.Machine$double.eps)
  piGrid = seq(piEps, 1 - piEps, length.out = nPiGrid)
  shapeGrid = seq(prior$range[1], prior$range[2], length.out = nShapeGrid)

  logPosterior = outer(
    piGrid,
    shapeGrid,
    Vectorize(function(pi, shape) {
      zizLogLikelihood(obsData, counts, pi, shape) +
        dbeta(pi, shape1, shape2, log = TRUE) +
        prior$logd(shape)
    })
  )

  finiteValues = logPosterior[is.finite(logPosterior)]
  if (length(finiteValues) == 0L) {
    stop("The zero-inflated zeta posterior grid has no finite posterior values")
  }

  # Subtract the maximum log posterior before exponentiating. This standard
  # log-sum-exp rescaling preserves relative weights while avoiding overflow
  # and severe underflow during numerical integration.
  logScale = max(finiteValues)
  scaledWeights = exp(logPosterior - logScale)
  scaledWeights[!is.finite(scaledWeights)] = 0

  # The grid is rectangular and equally spaced, so each cell has integration
  # weight dPi * dShape. Marginal densities below retain the complementary
  # spacing factor so subsequent one-dimensional sums approximate integrals.
  dPi = diff(range(piGrid)) / (length(piGrid) - 1L)
  dShape = diff(range(shapeGrid)) / (length(shapeGrid) - 1L)
  scaledIntegral = sum(scaledWeights) * dPi * dShape

  if (!is.finite(scaledIntegral) || scaledIntegral <= 0) {
    stop("The zero-inflated zeta posterior grid could not be normalized")
  }

  jointDensity = scaledWeights / scaledIntegral
  marginalPiDensity = rowSums(jointDensity) * dShape
  marginalShapeDensity = colSums(jointDensity) * dPi

  piMean = sum(piGrid * marginalPiDensity) * dPi
  shapeMean = sum(shapeGrid * marginalShapeDensity) * dShape
  piVariance = sum((piGrid - piMean)^2 * marginalPiDensity) * dPi
  shapeVariance = sum((shapeGrid - shapeMean)^2 * marginalShapeDensity) * dShape
  covariance = sum(
    outer(piGrid - piMean, shapeGrid - shapeMean) * jointDensity
  ) * dPi * dShape

  list(
    pi = piGrid,
    shape = shapeGrid,
    logPosterior = logPosterior,
    density = jointDensity,
    marginalDensity = list(
      pi = marginalPiDensity,
      shape = marginalShapeDensity
    ),
    normalizingConstant = exp(logScale) * scaledIntegral,
    dPi = dPi,
    dShape = dShape,
    mean = c(pi = piMean, shape = shapeMean),
    varCov = matrix(
      c(piVariance, covariance, covariance, shapeVariance),
      nrow = 2L,
      dimnames = list(c("pi", "shape"), c("pi", "shape"))
    )
  )
}

# ---- zizPosteriorImportance.R ----
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

# ---- zizPosteriorLaplace.R ----
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

# ---- zizPosteriorProbabilities.R ----
#' Summarise posterior zero-inflated zeta probabilities
#'
#' Converts posterior parameter representations into a common set of posterior
#' summaries for P or S probabilities.
#'
#' @param probabilities Numeric matrix with posterior realisations in rows and
#'   probability terms in columns.
#' @param weights Optional non-negative posterior weights.
#' @param level Equal-tailed credible interval level.
#' @param posteriorMethod Posterior engine used to produce the representation.
#'
#' @return A data frame containing posterior means, standard deviations, and
#'   equal-tailed credible intervals.
#'
#' @keywords internal
summariseZizProbabilities = function(probabilities,
                                     weights = NULL,
                                     level = 0.95,
                                     posteriorMethod = NA_character_) {
  probabilities = as.matrix(probabilities)

  if (!is.numeric(probabilities) || nrow(probabilities) == 0L ||
      ncol(probabilities) == 0L || any(!is.finite(probabilities))) {
    stop("probabilities must be a non-empty finite numeric matrix")
  }

  if (is.null(colnames(probabilities)) || any(colnames(probabilities) == "")) {
    stop("probabilities must have non-empty column names")
  }

  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("level must be one finite number strictly between zero and one")
  }

  if (is.null(weights)) {
    weights = rep(1 / nrow(probabilities), nrow(probabilities))
  } else {
    if (!is.numeric(weights) || length(weights) != nrow(probabilities) ||
        any(!is.finite(weights)) || any(weights < 0)) {
      stop("weights must be finite, non-negative, and match the number of rows")
    }

    weightSum = sum(weights)
    if (!is.finite(weightSum) || weightSum <= 0) {
      stop("weights must have a positive finite sum")
    }
    weights = weights / weightSum
  }

  lowerProbability = (1 - level) / 2
  upperProbability = 1 - lowerProbability
  estimates = colSums(probabilities * weights)
  secondMoments = colSums(probabilities^2 * weights)
  variances = pmax(0, secondMoments - estimates^2)

  intervals = vapply(
    seq_len(ncol(probabilities)),
    function(columnIndex) {
      weightedZizQuantile(
        values = probabilities[, columnIndex],
        weights = weights,
        probabilities = c(lowerProbability, upperProbability)
      )
    },
    numeric(2L)
  )

  data.frame(
    term = colnames(probabilities),
    estimate = unname(estimates),
    sd = sqrt(unname(variances)),
    lower = unname(intervals[1L, ]),
    upper = unname(intervals[2L, ]),
    level = rep(level, ncol(probabilities)),
    posteriorMethod = rep(posteriorMethod, ncol(probabilities)),
    stringsAsFactors = FALSE
  )
}

#' Weighted quantiles for posterior probability summaries
#'
#' @param values Numeric values.
#' @param weights Non-negative weights.
#' @param probabilities Quantile probabilities in `[0, 1]`.
#'
#' @return Numeric weighted quantiles.
#'
#' @keywords internal
weightedZizQuantile = function(values, weights, probabilities) {
  if (!is.numeric(values) || !is.numeric(weights) ||
      length(values) != length(weights) || length(values) == 0L ||
      any(!is.finite(values)) || any(!is.finite(weights)) || any(weights < 0)) {
    stop("values and weights must be finite numeric vectors of equal length")
  }

  if (!is.numeric(probabilities) || any(!is.finite(probabilities)) ||
      any(probabilities < 0 | probabilities > 1)) {
    stop("probabilities must lie in [0, 1]")
  }

  weightSum = sum(weights)
  if (!is.finite(weightSum) || weightSum <= 0) {
    stop("weights must have a positive finite sum")
  }

  orderIndex = order(values)
  sortedValues = values[orderIndex]
  sortedWeights = weights[orderIndex] / weightSum
  cumulativeWeights = cumsum(sortedWeights)

  vapply(
    probabilities,
    function(probability) {
      if (probability <= 0) {
        return(sortedValues[1L])
      }
      if (probability >= 1) {
        return(sortedValues[length(sortedValues)])
      }

      quantileIndex = which(cumulativeWeights >= probability)[1L]
      sortedValues[quantileIndex]
    },
    numeric(1L)
  )
}

#' Return the support indices used to label fitted P or S probabilities.
#'
#' @param type P- or S-survey type.
#' @param nterms Number of fitted P/S probability terms to retain.
#' @return Integer support indices for naming P/S probability terms.
#' @keywords internal
#' @noRd
zizProbabilityIndices = function(type, nterms) {
  type = match.arg(type, c("P", "S"))
  nterms = as.integer(nterms)

  if (!is.finite(nterms) || length(nterms) != 1L || nterms <= 0L) {
    stop("nterms must be a positive integer")
  }

  if (type == "P") {
    0:(nterms - 1L)
  } else {
    seq_len(nterms)
  }
}

#' Summarise fitted P/S probabilities over a numerical ZIZ posterior grid.
#'
#' @param posteriorGrid Numerical posterior-grid representation.
#' @param type P- or S-survey type.
#' @param nterms Number of fitted P/S probability terms to retain.
#' @param level Probability level for intervals or summaries.
#' @return A data frame of posterior fitted-probability summaries.
#' @keywords internal
#' @noRd
summariseZizGridProbabilities = function(posteriorGrid,
                                          type,
                                          nterms,
                                          level = 0.95) {
  piValues = rep(posteriorGrid$pi, times = length(posteriorGrid$shape))
  shapeValues = rep(posteriorGrid$shape, each = length(posteriorGrid$pi))
  weights = as.vector(posteriorGrid$density) *
    posteriorGrid$dPi * posteriorGrid$dShape

  probabilities = zizProbabilities(
    pi = piValues,
    shape = shapeValues,
    n = zizProbabilityIndices(type, nterms),
    type = type
  )

  summariseZizProbabilities(
    probabilities = probabilities,
    weights = weights,
    level = level,
    posteriorMethod = "numerical"
  )
}

#' Summarise fitted P/S probabilities over sampled ZIZ posterior representations.
#'
#' @param pi Zero-inflation probability on the natural scale.
#' @param shape Zeta shape parameter on the fitPS scale.
#' @param type P- or S-survey type.
#' @param nterms Number of fitted P/S probability terms to retain.
#' @param weights Non-negative normalised or normalisable weights.
#' @param level Probability level for intervals or summaries.
#' @param posteriorMethod Character label identifying the posterior approximation method.
#' @return A data frame of posterior fitted-probability summaries.
#' @keywords internal
#' @noRd
summariseZizSampleProbabilities = function(pi,
                                            shape,
                                            type,
                                            nterms,
                                            weights = NULL,
                                            level = 0.95,
                                            posteriorMethod) {
  probabilities = zizProbabilities(
    pi = pi,
    shape = shape,
    n = zizProbabilityIndices(type, nterms),
    type = type
  )

  summariseZizProbabilities(
    probabilities = probabilities,
    weights = weights,
    level = level,
    posteriorMethod = posteriorMethod
  )
}

# ---- zizProbabilities.R ----
#' Compute zero-inflated zeta P or S probabilities
#'
#' Maps one or more pairs of zero-inflated zeta parameters to the requested
#' P- or S-survey probabilities. This is the shared internal transformation
#' used by posterior probability summaries.
#'
#' @param pi Numeric vector of mixing probabilities in the interval `[0, 1]`.
#' @param shape Numeric vector of standard zeta shape parameters greater than
#'   one.
#' @param n Non-negative integer P-term indices or positive integer S-term
#'   indices.
#' @param type Survey type, either `"P"` or `"S"`.
#'
#' @return A numeric matrix with one row per parameter pair and one column per
#'   requested probability. Columns are named with the corresponding P or S
#'   term.
#'
#' @keywords internal
zizProbabilities = function(pi, shape, n, type) {
  if (!is.numeric(pi) || length(pi) == 0L || any(!is.finite(pi)) ||
      any(pi < 0 | pi > 1)) {
    stop("pi must contain finite values in [0, 1]")
  }

  if (!is.numeric(shape) || length(shape) == 0L || any(!is.finite(shape))) {
    stop("shape must contain finite numeric values")
  }
  validateZetaShape(shape)

  parameterCount = max(length(pi), length(shape))
  if (!length(pi) %in% c(1L, parameterCount) ||
      !length(shape) %in% c(1L, parameterCount)) {
    stop("pi and shape must have equal lengths or one must have length one")
  }

  pi = rep(pi, length.out = parameterCount)
  shape = rep(shape, length.out = parameterCount)

  type = normaliseSurveyType(type)
  n = normaliseProbabilityIndices(n, type)
  latentValues = latentPsValues(n, type)
  inflatedIndex = if (type == "P") 0L else 1L

  probabilities = vapply(
    seq_along(latentValues),
    function(termIndex) {
      values = (1 - pi) * dzetaStandard(
        latentValues[termIndex],
        shape = shape
      )
      if (n[termIndex] == inflatedIndex) {
        values = values + pi
      }
      values
    },
    numeric(parameterCount)
  )

  if (parameterCount == 1L) {
    probabilities = matrix(probabilities, nrow = 1L)
  }

  colnames(probabilities) = psProbabilityTermNames(n, type)
  probabilities
}

# ---- psModel.R:zizModel ----
#' @rdname zetaModel
#' @export
zizModel = function() {
  newPsModel(
    model = "ziz",
    parameterNames = c("pi", "shape"),
    supportedEngines = c("numerical", "mcmc", "laplace", "importance"),
    subclass = "zizModel"
  )
}

# ---- psModel.R:modelProbabilities.zizModel ----
#' @rdname modelProbabilities
#' @keywords internal
#' @exportS3Method modelProbabilities zizModel
#' @noRd
modelProbabilities.zizModel = function(model, parameters, n, type, ...) {
  zizProbabilities(
    pi = modelParameter(parameters, "pi"),
    shape = modelParameter(parameters, "shape"),
    n = n,
    type = type
  )
}

# ---- psModel.R:modelLogLikelihood.zizModel ----
#' Evaluate the zero-inflated zeta log likelihood.
#'
#' @param model A `zizModel` descriptor.
#' @param parameters Named parameters containing scalar `pi` and `shape`.
#' @param data An object of class `psData`.
#' @param ... Additional arguments reserved for future ZIZ likelihood controls.
#' @return Scalar log likelihood.
#' @keywords internal
#' @exportS3Method modelLogLikelihood zizModel
#' @noRd
modelLogLikelihood.zizModel = function(model, parameters, data, ...) {
  if (!is(data, "psData")) {
    stop("data must be an object of class psData")
  }

  pi = modelParameter(parameters, "pi")
  shape = modelParameter(parameters, "shape")

  if (!is.numeric(pi) || length(pi) != 1L || !is.finite(pi)) {
    stop("pi must be one finite numeric value")
  }
  if (!is.numeric(shape) || length(shape) != 1L || !is.finite(shape)) {
    stop("shape must be one finite numeric value")
  }

  observations = modelObservationData(model, data)
  zizLogLikelihood(
    obsData = observations,
    counts = data$data$rn,
    pi = pi,
    shape = shape
  )
}

# ---- fit.R:fitModel.zizModel ----
#' @rdname fitModel
#' @keywords internal
#' @exportS3Method fitModel zizModel
#' @noRd
fitModel.zizModel = function(model, x, ...) {
  result = fitZIDistImpl(x = x, ...)
  result$modelObject = model
  result
}

# ---- zizLaplacePosteriorEngine.R:fitLaplacePosteriorModel.zizModel ----
#' Fit a Laplace posterior for the zero-inflated zeta model.
#'
#' This method preserves the established transformed-coordinate optimisation,
#' Hessian calculation, and delta-method covariance while wrapping the result
#' in the common posterior-engine representation.
#'
#' @param model A `zizModel` descriptor.
#' @param engine A `laplacePosteriorEngine` object.
#' @param x An object of class `psData`.
#' @param prior A prior object created by `makePrior()`.
#' @param shape1 First beta-prior shape parameter for `pi`.
#' @param shape2 Second beta-prior shape parameter for `pi`.
#' @param start Starting values for `pi` and `shape`.
#' @param ... Additional Laplace controls reserved for future use.
#' @return A `laplacePosteriorRepresentation` object.
#' @keywords internal
#' @exportS3Method fitLaplacePosteriorModel zizModel
#' @noRd
fitLaplacePosteriorModel.zizModel = function(model,
                                              engine,
                                              x,
                                              prior,
                                              shape1 = 1,
                                              shape2 = 1,
                                              start = c(pi = 0.5, shape = 2),
                                              ...) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData")
  }
  if (!inherits(prior, "psPrior")) {
    stop("prior must be an object of class psPrior")
  }

  validateZetaPriorRange(prior$range)

  modelObservationData(model, x)
  approximation = makeZizPosteriorLaplace(
    x = x,
    prior = prior,
    shape1 = shape1,
    shape2 = shape2,
    start = start
  )

  newPsPosteriorRepresentation(
    engine = engine,
    value = list(
      approximation = approximation,
      mean = approximation$mode,
      variance = approximation$varCov
    ),
    metadata = list(
      model = model$model,
      convergence = approximation$convergence,
      logPosteriorMode = approximation$logPosteriorMode,
      shape1 = shape1,
      shape2 = shape2
    )
  )
}

# ---- zizImportancePosteriorEngine.R:fitImportancePosteriorModel.zizModel ----
#' Fit an importance-sampling posterior for the zero-inflated zeta model.
#'
#' This method preserves the established Laplace-based Gaussian proposal,
#' transformed-coordinate sampling, and normalized importance weights while
#' wrapping the result in the common posterior-engine representation.
#'
#' @param model A `zizModel` descriptor.
#' @param engine An `importancePosteriorEngine` object.
#' @param x An object of class `psData`.
#' @param prior A prior object created by `makePrior()`.
#' @param shape1 First beta-prior shape parameter for `pi`.
#' @param shape2 Second beta-prior shape parameter for `pi`.
#' @param nSamples Number of importance samples.
#' @param proposalScale Positive proposal covariance scale multiplier.
#' @param seed Optional random-number seed.
#' @param start Starting values for `pi` and `shape`.
#' @param laplace Optional precomputed Laplace approximation.
#' @param ... Additional importance controls reserved for future use.
#' @return An `importancePosteriorRepresentation` object.
#' @keywords internal
#' @exportS3Method fitImportancePosteriorModel zizModel
#' @noRd
fitImportancePosteriorModel.zizModel = function(model,
                                                 engine,
                                                 x,
                                                 prior,
                                                 shape1 = 1,
                                                 shape2 = 1,
                                                 nSamples = 5000,
                                                 proposalScale = 2,
                                                 seed = NULL,
                                                 start = c(pi = 0.5, shape = 2),
                                                 laplace = NULL,
                                                 ...) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData")
  }
  if (!inherits(prior, "psPrior")) {
    stop("prior must be an object of class psPrior")
  }

  validateZetaPriorRange(prior$range)

  modelObservationData(model, x)
  approximation = makeZizPosteriorImportance(
    x = x,
    prior = prior,
    shape1 = shape1,
    shape2 = shape2,
    nSamples = nSamples,
    proposalScale = proposalScale,
    seed = seed,
    start = start,
    laplace = laplace
  )

  newPsPosteriorRepresentation(
    engine = engine,
    value = list(
      approximation = approximation,
      mean = approximation$mean,
      variance = approximation$varCov
    ),
    metadata = c(
      list(
        model = model$model,
        shape1 = shape1,
        shape2 = shape2,
        seed = seed
      ),
      approximation$diagnostics
    )
  )
}

# ---- bayesianPosteriorOrchestration.R:summariseNumericalPosteriorProbabilities.zizModel ----
#' @rdname summariseNumericalPosteriorProbabilities
#' @keywords internal
#' @exportS3Method summariseNumericalPosteriorProbabilities zizModel
#' @noRd
summariseNumericalPosteriorProbabilities.zizModel = function(model,
                                                              engine,
                                                              representation,
                                                              x,
                                                              nterms,
                                                              level = 0.95,
                                                              ...) {
  summariseZizGridProbabilities(
    posteriorGrid = representation$value$posteriorGrid,
    type = x$type,
    nterms = nterms,
    level = level
  )
}
