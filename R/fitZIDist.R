#' Fit a Zero-Inflated Zeta Distribution to Forensic Data
#'
#' This function uses maximum likelihood estimation (MLE) or Bayesian estimation
#' (MCMC) to estimate the mixing parameter and the shape parameter of a
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
#' @details The function returns an object of class \code{psFit} which is a
#'   \code{list} containing eight or nine elements:
#' \describe{
#' \item{\code{psData}}{ -- an object of class \code{psData}--see \code{\link{readData}},}
#' \item{\code{fit}}{ -- the fitted object from \code{\link[stats]{optim}},}
#' \item{\code{pi}}{ - the maximum likelihood estimate, or the posterior mean, of the mixing parameter,}
#' \item{\code{shape}}{ -- the maximum likelihood estimate, or the posterior mean, of the shape parameter,}
#' \item{\code{var.cov}}{ -- the estimated (posterior) variance-covariance matrix for the parameters,}
#' \item{\code{fitted}}{ -- a named \code{vector} containing the first \code{nterms} of the fitted distribution.}
#' \item{\code{model}}{ -- set to \code{"ziz"} for this model,}
#' \item{\code{method}}{ -- the method of estimation used, either \code{"mle"} or \code{"bayes"},}
#' \item{\code{chain}}{ -- if \code{method == "bayes"}, then this element will contain the Markov Chain from the sampler,
#' that is, hopefully a sample from the posterior density of the mixing parameter and the shape parameter.
#' If \code{method == "mle"}, then this element does not exist.}
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
#' \code{makePrior()} is used.

#'
#' @seealso \code{\link{plot.psFit}}, \code{\link{print.psFit}},
#'   \code{\link{probfun}}.
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
#'   If \code{method == "bayes"}, then there are seven optional parameters (which,
#'   despite the documentation, are actually case-insensitive):
#' \describe{
#'   \item{\code{theta0}}{ -- The initial values of the mixing parameter and shape parameter. The default is \code{c(0.5, 2)}. }
#'   \item{\code{a}}{ -- The lower bound for the default uniform prior on \eqn{\log(\mathrm{shape} - 1)}{log(shape - 1)}. The default is -2. }
#'   \item{\code{b}}{ -- The upper bound for the default uniform prior on \eqn{\log(\mathrm{shape} - 1)}{log(shape - 1)}. The default is +2. }
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
#' @return an object of class \code{psFit}--see Details.
#' @export
#'
#' @examples
#' data(Psurveys)
#' roux = Psurveys$roux
#' fit = fitZIDist(roux)
#' fit
fitZIDist = function(x, nterms = 10,
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
      model = "ziz",
      method = "mle"
    )


    class(result) = "psFit"

    return(result)
  }else{ ## method == "bayes"
    options = if (missing(prior)) {
      normaliseBayesOptions(bayesOptions = bayesOptions)
    } else {
      normaliseBayesOptions(bayesOptions = bayesOptions, prior = prior)
    }

    if (options$posteriorMethod == "numerical") {
      result = fitZIDistBayesNumerical(
        x = x,
        nterms = nterms,
        prior = options$prior,
        ...
      )
      result$bayesOptions = options
      return(result)
    }

    if (options$posteriorMethod == "mcmc") {
      result = fitZIDistBayes(x = x, nterms = nterms, ...)
      result$posteriorMethod = "mcmc"
      result$bayesOptions = options
      return(result)
    }

    if (options$posteriorMethod == "laplace") {
      result = fitZIDistBayesLaplace(
        x = x,
        nterms = nterms,
        prior = options$prior,
        ...
      )
      result$bayesOptions = options
      return(result)
    }

    if (options$posteriorMethod == "importance") {
      result = fitZIDistBayesImportance(
        x = x,
        nterms = nterms,
        prior = options$prior,
        ...
      )
      result$bayesOptions = options
      return(result)
    }

    stop(
      "posteriorMethod = ",
      sQuote(options$posteriorMethod),
      " is not implemented for fitZIDist() yet"
    )
  }
}

#' @rdname fitZIDist
#' @export
fitZIdist = fitZIDist

#' @rdname fitZIDist
#' @export
fitzidist = fitZIDist

zi.loglik = function(y, theta){
  p = theta[1]
  shape = theta[2]
  rval = (1 - p) * dzetaStandard(y$n, shape = shape)
  rval[y$n == 1] = rval[y$n == 1] + p
  sum(y$rn *log(rval))
}


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
    model = "ziz"
  )


  class(result) = "psFit"

  return(result)
}

