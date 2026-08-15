# Zeta model implementation
#
# Distribution-specific zeta fitting, probability, posterior, simulation, and
# model-dispatch code is intentionally consolidated in this file.

# ---- fitDist.R ----
#' Fit a Zeta Distribution to Forensic Data
#'
#' This function uses maximum likelihood estimation (MLE), or Bayesian
#' estimation (MCMC), to estimate the shape parameter of a zeta distribution
#' from a set of observed counts for either the number of groups/sources of
#' forensically interesting material (mostly glass or paint) recovered from
#' clothing, or the number of fragments/particles in each group. This, in turn,
#' allows the estimation of the P and S probabilities, as described by Evett and
#' Buckleton (1990), which are used in computing the likelihood ratio (LR) for
#' activity level propositions. The data arise from clothing surveys.
#' The general method is described in Coulson et al. (2001), although poor
#' typesetting and a lack of defined terms make it hard to follow. This
#' package improves on the estimation in that linear interpolation is not
#' required, and standard numerical optimisation is used instead. The zeta
#' distribution has probability mass function \deqn{p(k) =
#' \frac{k^{-s}}{\zeta(s)}}{p(k) = k^-s/zeta(s)} where \eqn{\zeta(s)}{zeta(s)}
#' is the Riemann zeta function. Coulson et al. (2001) did not have an easy way
#' to rapidly compute this quantity, hence their use of linear interpolation.
#'
#' @aliases fitdist
#'
#' @details The function returns an object of class \code{psFit}. Core fitted
#'   values are stored at the top level, while Bayesian uncertainty is stored
#'   in the attached \code{psPosterior} object. Important components include:
#' \describe{
#' \item{\code{psData}}{ -- an object of class \code{psData}--see \code{\link{readData}},}
#' \item{\code{fit}}{ -- method-specific point-estimate information,}
#' \item{\code{shape}}{ -- the maximum likelihood estimate, or the posterior mean, of the shape parameter,}
#' \item{\code{var.shape}}{ -- the maximum likelihood estimate, or posterior estimate, of the variance of the shape parameter,}
#' \item{\code{fitted}}{ -- a named \code{vector} containing the first \code{nterms} of
#' the fitted distribution.}
#' \item{\code{model}}{ -- set to \code{"zeta"} for this model.}
#' \item{\code{method}}{ -- the method of estimation used, either \code{"mle"} or \code{"bayes"}.}
#' \item{\code{posterior}}{ -- for Bayesian fits, a \code{psPosterior} object containing parameter summaries, probability summaries, diagnostics, and the engine-specific representation,}
#' \item{\code{chain}}{ -- if \code{method == "bayes"}, then this element will contain the Markov Chain from the sampler,
#' that is, hopefully a sample from the posterior density of the shape parameter. If \code{method == "mle"}, then this element does not exist.}
#' }
#'   The output can be used in a variety of ways. If the interest is just in the
#'   shape parameter estimate, then the \code{shape} member of the \code{psFit}
#'   object contains this information. It is also displayed along with a number
#'   of fitted probabilities by the \code{\link{print.psFit}} method. The fitted
#'   object can also be plotted using the plot method \code{\link{plot.psFit}},
#'   and to create a probability function with \code{\link{probfun}}. The
#'   \code{shape} value stored in the fitted object is the zeta distribution
#'   shape parameter and must satisfy \code{shape > 1}.
#'
#'   This function implements both maximum likelihood estimation (MLE) and
#'   Bayesian estimation. Both modes of estimation require additional information
#'   such as starting values and parameters for priors. Please read the
#'   documentation for the \code{...} argument closely because it explains what
#'   you can change and what the default values are.
#'
#'   Bayesian estimation uses the prior returned by \code{\link{makePrior}}
#'   and the posterior engine selected by \code{bayesOptions$posteriorMethod}.
#'   Plain zeta currently supports \code{"numerical"} and \code{"mcmc"}.
#'   Both engines return the same \code{psPosterior} contract, including
#'   posterior parameter summaries and posterior summaries of the fitted P or S
#'   probabilities.
#'
#' @seealso \code{\link{plot.psFit}}, \code{\link{print.psFit}},
#'   \code{\link{probfun}}.
#'
#' @references Coulson, S. A., Buckleton, J. S., Gummer, A. B., and Triggs,
#'   C.M., "Glass on clothing and shoes of members of the general population and
#'   people suspected of breaking crimes", Science & Justice 2001: 41(1):
#'   39--48.
#'
#'   Evett, I. W. and Buckleton, J. S., "The interpretation of glass evidence. A
#'   practical approach", Journal of the Forensic Science Society 1990: 30(4):
#'   215--223.
#'
#' @param x an object of type \code{psData}, usually obtained from
#'   \code{\link{readData}}.
#' @param nterms the number of terms to compute the probability distribution
#'   for.
#' @param method primary fitting method. Use \code{"mle"} for maximum
#'   likelihood estimation or \code{"bayes"} for Bayesian estimation. Legacy
#'   Bayesian aliases \code{"integrate"}, \code{"numerical"}, and
#'   \code{"mcmc"} are accepted with a deprecation warning and translated to
#'   \code{method = "bayes"} with the corresponding
#'   \code{bayesOptions$posteriorMethod}.
#' @param prior optional prior object used by the Bayesian methods. This is
#'   retained for backward compatibility. New code should usually pass priors
#'   through \code{bayesOptions}. If omitted, \code{makePrior()} is used.
#' @param bayesOptions optional list controlling Bayesian fitting. The
#'   \code{posteriorMethod} element selects \code{"numerical"},
#'   \code{"mcmc"}, \code{"laplace"}, or \code{"importance"}. The
#'   default is \code{"numerical"}. The \code{prior} element may contain
#'   a prior object returned by \code{makePrior()}.
#' @param ... other arguments that control the estimation methods. If
#'   \code{method == "mle"}, then the user can provide an optional argument
#'   \code{start} which is the starting value for the numerical optimisation. If
#'   this is not provided, then \code{start = 1} by default. If you specify your
#'   own starting value, it must satisfy \code{shape > 1}.
#'
#'   If \code{method == "bayes"}, then there are five optional parameters (which,
#'   despite the documentation, are actually case-insensitive):
#' \describe{
#'  \item{\code{shape0}}{ -- The initial value of the shape parameter. The default is 2.}
#'  \item{\code{a}}{ -- The lower bound for the default uniform prior on \eqn{\log(\mathrm{shape} - 1)}{log(shape - 1)}. The default is -2.}
#'  \item{\code{b}}{ -- The upper bound for the default uniform prior on \eqn{\log(\mathrm{shape} - 1)}{log(shape - 1)}. The default is +2.}
#'  \item{\code{nIter}}{ -- The number of samples to save from the chain. Must be greater than zero, and ideally greater than 1000.}
#'  \item{\code{nBurnIn}}{ -- The number of samples to discard from the chain. Must be greater than zero. **NOTE**: the sampler runs for \code{nIter + nBurnIn} iterations,
#'  so you do not need to factor this number into your number of samples, \code{nIter}.}
#'  \item{\code{silent}}{ -- A logical variable which allows the user to get a progress bar if they want. \code{TRUE} by default.}
#' }
#'
#' @importFrom stats density integrate optim runif splinefun
#' @importFrom VGAM dzeta
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @section Deprecated:
#' `fitDist()` is retained as a compatibility wrapper. New code should use
#' `fit(x, model = zetaModel())`.
#'
#' @return an object of class \code{psFit}--see Details.
#' @export
#'
#' @examples
#' p = readData(system.file("extdata", "p.xlsx", package = "fitPS"))
#' fit = fitDist(p)
#' fit
#'
#' ## Compare to the Bayesian estimates
#' fit2 = fitDist(p, method = "bayes")
#' fit2
#'
#' fit3 = fitDist(
#'   p,
#'   method = "bayes",
#'   bayesOptions = list(posteriorMethod = "numerical")
#' )
#' fit3
fitDist = function(x,
                    nterms = 10,
                    method = c(
                      "mle", "bayes", "integrate", "numerical", "mcmc",
                      "laplace", "importance"
                    ),
                    prior,
                    bayesOptions = NULL,
                    ...) {
  signalDeprecatedFitter(
    old = "fitDist",
    replacement = "fit(x, model = zetaModel())"
  )

  methodInfo = normaliseBayesMethod(method, bayesOptions = bayesOptions)
  args = list(
    x = x,
    model = zetaModel(),
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

#' Internal zeta fitting implementation.
#'
#' @inheritParams fitDist
#' @return An object of class `psFit`.
#' @keywords internal
#' @noRd
fitDistImpl = function(x, nterms = 10,
                   method = c("mle", "bayes", "integrate", "numerical", "mcmc", "laplace", "importance"),
                   prior,
                   bayesOptions = NULL,
                   ...){
  nvals = 1:nterms
  if(!is(x, "psData")){
    stop("x must be an object of class psData")
  }

  obsData = if(x$type == 'P'){ ## the main difference is that the values need 1 added
                x$data$n + 1
            }else{
              x$data$n
            }

  methodInfo = normaliseBayesMethod(method, bayesOptions = bayesOptions)
  method = methodInfo$method
  bayesOptions = methodInfo$bayesOptions

  model = zetaModel()

  if(method == "mle"){
    validateMleObservationSupport(x)

    dotargs = list(...)
    if("start" %in% names(dotargs)){
      start = dotargs$start
    }else if("shape" %in% names(dotargs)){
      start = dotargs$shape
    }else{
      start = 2
    }

    validateZetaShape(start, "start")

    logLik = function(shape){
      -sum(x$data$rn * dzetaStandard(obsData, shape = shape, log = TRUE))
    }

    # fit = nlminb(start = start,
    #              objective = logLik,
    #              lower = 1)

    fit = optim(par = start,
                  fn = logLik,
                  method = "L-BFGS-B",
                  lower = 1 + sqrt(.Machine$double.eps),
                  hessian = TRUE)

    shape = fit$par
    N = sum(x$data$rn)

    vshape = function(shape){
      z = VGAM::zeta(shape)
      zprime = VGAM::zeta(shape, 1)
      zprimeprime = VGAM::zeta(shape, 2)
      numer = z^2
      denom = N *(z * zprimeprime - zprime^2)

      return(numer / denom)
    }

    var.shape = vshape(shape)

    fitted = dzetaStandard(nvals, shape = shape)
    names(fitted) = if(x$type == 'P'){
      paste0("P", nvals - 1)
    }else{
      paste0("S", nvals)
    }

    result = list(
      psData = x,
      fit = fit,
      shape = shape,
      var.shape = var.shape,
      fitted = fitted,
      model = model$model,
      modelObject = model,
      method = "mle"
    )

    class(result) = "psFit"

    return(result)
  } else if (method == "bayes") {
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
  } else {
    stop("Unknown method:", method)
  }
}

#' @describeIn fitDist Fit a Zeta Distribution to Forensic Data
#' @export
fitdist = fitDist

# ---- internalfunctions.R ----
#' Evaluate the zeta log likelihood for internal Bayesian calculations.
#'
#' @param x An input object or numeric vector required by the helper.
#' @param shape Zeta shape parameter on the fitPS scale.
#' @return A numeric log-likelihood value.
#' @keywords internal
#' @noRd
zetaloglikelihood = function(x, shape){
  offset = ifelse(x$type == 'P', 1, 0)
  if(length(shape) > 1){
    ll = function(s){
      sum(x$data$rn * dzetaStandard(x$data$n + offset, shape = s, log = TRUE))
    }
    sapply(shape, ll)
  }else{
    sum(x$data$rn * dzetaStandard(x$data$n + offset, shape = shape, log = TRUE))
  }
}

zll = zetaloglikelihood

#' Evaluate the zeta likelihood on the natural probability scale.
#'
#' @param x An input object or numeric vector required by the helper.
#' @param shape Zeta shape parameter on the fitPS scale.
#' @return A numeric likelihood value.
#' @keywords internal
#' @noRd
zetalikelihood = function(x, shape){
  exp(zll(x, shape))
}

#' Evaluate the legacy bounded uniform prior density for the zeta shape.
#'
#' @param s Candidate zeta shape value.
#' @param a Lower bound of the legacy uniform prior.
#' @param b Upper bound of the legacy uniform prior.
#' @return A numeric prior density.
#' @keywords internal
#' @noRd
zetaunifprior = function(s, a = exp(-2), b = exp(2)){
  result = numeric(length(s))
  inRange = s >= a & s <= b
  result[inRange] = 1 / ((b - a) * s[inRange])
  result[!inRange] = 0
  return(result)
}

zup = zetaunifprior

# ---- rzeta.R ----
#' Generate random variates from a zeta distribution
#'
#' @param n Same as \code{\link[stats]{Poisson}}.
#' @param shape The standard zeta shape parameter, greater than 1.
#'
#' See \code{\link[VGAM]{rzeta}}.
#'
#' @export
rzeta = function(n, shape){
  rzetaStandard(n = n, shape = shape)
}

# ---- zetaParameterisation.R ----
#' Validate that a zeta shape parameter is finite and greater than one.
#'
#' @param shape Zeta shape parameter on the fitPS scale.
#' @param argument Argument name used in validation error messages.
#' @return The validated shape value, invisibly or unchanged as used by callers.
#' @keywords internal
#' @noRd
validateZetaShape = function(shape, argument = "shape"){
  if(any(!is.finite(shape))){
    stop(argument, " must be finite.")
  }

  if(any(shape <= 1)){
    stop(argument, " must be greater than 1.")
  }

  invisible(shape)
}

#' Convert the fitPS zeta shape convention to the VGAM parameter convention.
#'
#' @param shape Zeta shape parameter on the fitPS scale.
#' @return The equivalent VGAM zeta parameter.
#' @keywords internal
#' @noRd
standardToVgamShape = function(shape){
  validateZetaShape(shape)
  shape - 1
}

#' Evaluate the zeta mass function using the fitPS shape convention.
#'
#' @param x An input object or numeric vector required by the helper.
#' @param shape Zeta shape parameter on the fitPS scale.
#' @param log Logical; return log probabilities when `TRUE`.
#' @return Numeric zeta probabilities or log probabilities.
#' @keywords internal
#' @noRd
dzetaStandard = function(x, shape, log = FALSE){
  vgamShape = standardToVgamShape(shape)
  VGAM::dzeta(x, shape = vgamShape, log = log)
}

#' Generate zeta random variates using the fitPS shape convention.
#'
#' @param n Number of random variates or observations.
#' @param shape Zeta shape parameter on the fitPS scale.
#' @return Integer-valued zeta random variates.
#' @keywords internal
#' @noRd
rzetaStandard = function(n, shape){
  vgamShape = standardToVgamShape(shape)
  VGAM::rzeta(n = n, shape = vgamShape)
}

# ---- zetaProbabilities.R ----
#' Compute zeta P or S probabilities
#'
#' Maps one or more standard zeta shape parameters to the requested P- or
#' S-survey probabilities. This is the frequentist bootstrap counterpart to
#' the zero-inflated transformation in `zizProbabilities()`.
#'
#' @param shape Numeric vector of standard zeta shape parameters greater than
#'   one.
#' @param n Non-negative integer P-term indices or positive integer S-term
#'   indices.
#' @param type Survey type, either `"P"` or `"S"`.
#'
#' @return A numeric matrix with one row per shape value and one column per
#'   requested probability. Columns are named with the corresponding P or S
#'   term.
#'
#' @keywords internal
zetaProbabilities = function(shape, n, type) {
  if (!is.numeric(shape) || length(shape) == 0L || any(!is.finite(shape))) {
    stop("shape must contain finite numeric values")
  }
  validateZetaShape(shape)

  type = normaliseSurveyType(type)
  n = normaliseProbabilityIndices(n, type)
  latentValues = latentPsValues(n, type)

  probabilities = vapply(
    latentValues,
    function(value) {
      dzetaStandard(value, shape = shape)
    },
    numeric(length(shape))
  )

  if (length(shape) == 1L) {
    probabilities = matrix(probabilities, nrow = 1L)
  }

  colnames(probabilities) = psProbabilityTermNames(n, type)
  probabilities
}

# ---- psModel.R:zetaModel ----
#' Construct built-in fitPS model descriptors
#'
#' These constructors create model objects for the built-in fitPS distributions.
#' They can be supplied directly to [fit()]. The lower-level third-party model
#' constructor and extension generics are introduced separately after the public
#' contract has been validated.
#'
#' @return An object inheriting from `psModel` and the concrete built-in model
#'   class.
#' @export
zetaModel = function() {
  newPsModel(
    model = "zeta",
    parameterNames = "shape",
    supportedEngines = c("numerical", "mcmc"),
    subclass = "zetaModel"
  )
}

# ---- psModel.R:modelProbabilities.zetaModel ----
#' @rdname modelProbabilities
#' @keywords internal
#' @exportS3Method modelProbabilities zetaModel
#' @noRd
modelProbabilities.zetaModel = function(model, parameters, n, type, ...) {
  zetaProbabilities(
    shape = modelParameter(parameters, "shape"),
    n = n,
    type = type
  )
}

# ---- psModel.R:modelLogLikelihood.zetaModel ----
#' Evaluate the plain zeta log likelihood.
#'
#' @param model A `zetaModel` descriptor.
#' @param parameters Named parameters containing scalar `shape`.
#' @param data An object of class `psData`.
#' @param ... Additional arguments reserved for future zeta likelihood controls.
#' @return Scalar log likelihood.
#' @keywords internal
#' @exportS3Method modelLogLikelihood zetaModel
#' @noRd
modelLogLikelihood.zetaModel = function(model, parameters, data, ...) {
  if (!is(data, "psData")) {
    stop("data must be an object of class psData")
  }

  shape = modelParameter(parameters, "shape")
  if (!is.numeric(shape) || length(shape) != 1L || !is.finite(shape)) {
    stop("shape must be one finite numeric value")
  }
  validateZetaShape(shape)

  observations = modelObservationData(model, data)
  sum(data$data$rn * dzetaStandard(observations, shape = shape, log = TRUE))
}

# ---- fit.R:fitModel.zetaModel ----
#' @rdname fitModel
#' @keywords internal
#' @exportS3Method fitModel zetaModel
#' @noRd
fitModel.zetaModel = function(model, x, ...) {
  result = fitDistImpl(x = x, ...)
  result$modelObject = model
  result
}

# ---- zetaMcmcPosterior.R:fitMcmcPosteriorModel.zetaModel ----
#' @rdname fitMcmcPosteriorModel
#' @keywords internal
#' @exportS3Method fitMcmcPosteriorModel zetaModel
#' @noRd
fitMcmcPosteriorModel.zetaModel = function(model,
                                            engine,
                                            x,
                                            prior,
                                            shape0 = 2,
                                            nIter = 1e4,
                                            nBurnIn = 1e3,
                                            silent = TRUE,
                                            ...) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData")
  }
  if (!inherits(prior, "psPrior")) {
    stop("prior must be an object of class psPrior")
  }

  validateZetaShape(shape0, "shape0")
  validateZetaPriorRange(prior$range)
  modelObservationData(model, x)

  if (nIter < 1000) {
    warning(
      "The number of samples from the MCMC chain really should be 1000 or higher."
    )
  }
  if (nIter <= 0 || nBurnIn <= 0) {
    stop("nIter and nBurnIn must be greater than zero.")
  }

  lower = prior$range[1]
  upper = prior$range[2]
  nTotal = nIter + nBurnIn

  # Draw the complete proposal and acceptance streams before the loop. Keeping
  # this ordering preserves the legacy RNG sequence for a fixed set.seed().
  proposals = runif(nTotal, lower, upper)
  logUniforms = log(runif(nTotal))
  chain = numeric(nIter)

  logPosterior = function(shape) {
    modelLogLikelihood(
      model = model,
      parameters = list(shape = shape),
      data = x
    ) + prior$logd(shape)
  }

  currentLogPosterior = logPosterior(shape0)
  if (!is.finite(currentLogPosterior)) {
    stop("Log likelihood is not finite at starting value")
  }

  if (!silent) {
    progress = txtProgressBar(
      min = 1,
      max = nTotal,
      initial = 1,
      style = 3,
      label = "Burning in"
    )
  }

  i = 1
  while (i <= nTotal) {
    proposedShape = proposals[i]
    proposedLogPosterior = logPosterior(proposedShape)

    if (proposedLogPosterior > currentLogPosterior ||
        logUniforms[i] < (proposedLogPosterior - currentLogPosterior)) {
      shape0 = proposedShape
      currentLogPosterior = proposedLogPosterior
    }

    if (i > nBurnIn) {
      chain[i - nBurnIn] = shape0
    }

    i = i + 1
    if (!silent) {
      if (i <= nBurnIn) {
        setTxtProgressBar(progress, i)
      } else if (i <= nTotal) {
        setTxtProgressBar(progress, i, label = "Sampling")
      }
    }
  }

  if (!silent) {
    close(progress)
  }

  posteriorMean = mean(chain)
  posteriorVariance = var(chain)

  newPsPosteriorRepresentation(
    engine = engine,
    value = list(
      chain = chain,
      mean = c(shape = posteriorMean),
      variance = matrix(
        posteriorVariance,
        nrow = 1L,
        dimnames = list("shape", "shape")
      )
    ),
    metadata = list(
      model = model$model,
      bounds = c(lower = lower, upper = upper),
      nIter = nIter,
      nBurnIn = nBurnIn,
      acceptance = mean(diff(chain) != 0)
    )
  )
}

# ---- bayesianPosteriorOrchestration.R:establishedBayesianFitFields.zetaModel ----
#' @rdname establishedBayesianFitFields
#' @keywords internal
#' @exportS3Method establishedBayesianFitFields zetaModel
#' @noRd
establishedBayesianFitFields.zetaModel = function(model, representation) {
  validatePosteriorRepresentation(representation)

  variance = unname(representation$value$variance["shape", "shape"])
  result = list(var.shape = variance)

  if (!is.null(representation$value$chain)) {
    chain = representation$value$chain
    densityEstimate = density(
      chain,
      from = representation$metadata$bounds[["lower"]]
    )
    result$chain = chain
    result$pdf = splinefun(densityEstimate$x, densityEstimate$y)
  } else if (is.function(representation$value$density)) {
    result$pdf = representation$value$density
  }

  result
}

# ---- Zeta prior validation ----
#' Validate prior support for a zeta shape parameter
#'
#' @param range Two-element prior support.
#' @return `range`, invisibly, when valid.
#' @keywords internal
#' @noRd
validateZetaPriorRange = function(range) {
  validatePriorRange(range)
  if (range[1L] <= 1) {
    stop("zeta prior range must have lower bound greater than 1")
  }
  invisible(range)
}
