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
#' @details The function returns an object of class \code{psFit} which is a
#'   \code{list} containing seven or eight elements:
#' \describe{
#' \item{\code{psData}}{ -- an object of class \code{psData}--see \code{\link{readData}},}
#' \item{\code{fit}}{ -- the fitted object from \code{\link[stats]{optim}},}
#' \item{\code{shape}}{ -- the maximum likelihood estimate, or the posterior mean, of the shape parameter,}
#' \item{\code{var.shape}}{ -- the maximum likelihood estimate, or posterior estimate, of the variance of the shape parameter,}
#' \item{\code{fitted}}{ -- a named \code{vector} containing the first \code{nterms} of
#' the fitted distribution.}
#' \item{\code{model}}{ -- set to \code{"zeta"} for this model.}
#' \item{\code{method}}{ -- the method of estimation used, either \code{"mle"} or \code{"bayes"}.}
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
#'   Currently the Bayesian estimation is done using the prior returned by
#'   \code{\link{makePrior}}. By default this is a Uniform[a, b] prior on
#'   \eqn{\log(\mathrm{shape} - 1)}{log(shape - 1)}, so the prior support always has \code{shape > 1}. This may become more
#'   flexible in the future. Similarly, the estimation is done using a simple
#'   Metropolis-Hastings sampler. It might be more efficient to sample through
#'   adaptive rejection sampling, but it is unclear whether it is worth the
#'   effort.
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
fitDist = function(x, nterms = 10,
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

  if(method == "mle"){

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
      model = "zeta",
      method = "mle"
    )

    class(result) = "psFit"

    return(result)
  }else if (method == "bayes") {
    options = if (missing(prior)) {
      normaliseBayesOptions(bayesOptions = bayesOptions)
    } else {
      normaliseBayesOptions(bayesOptions = bayesOptions, prior = prior)
    }

    if (options$posteriorMethod == "numerical") {
      result = fitDistBayesIntegrate(
        x = x,
        prior = options$prior,
        nterms = nterms,
        ...
      )
    } else if (options$posteriorMethod == "mcmc") {
      result = fitDistBayes(
        x = x,
        prior = options$prior,
        nterms = nterms,
        ...
      )
    } else {
      stop(
        "posteriorMethod = ",
        sQuote(options$posteriorMethod),
        " is not implemented for fitDist() yet"
      )
    }

    result$method = "bayes"
    result$posteriorMethod = options$posteriorMethod
    result$bayesOptions = options
    return(result)
  } else {
    stop("Unknown method:", method)
  }
}

#' @describeIn fitDist Fit a Zeta Distribution to Forensic Data
#' @export
fitdist = fitDist
