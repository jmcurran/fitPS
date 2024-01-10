#' Fit a Zero-Inflated Zeta Distribution to Forensic Data
#'
#' This function uses maximum likelihood estimation (MLE) to estimate mixing
#' parameter and the shape parameter of a zero-inflated zeta distribution from a
#' set of observed counts for either the number of groups/sources of
#' forensically interesting material (mostly glass or paint) recovered from
#' clothing, or the number of fragments/particles in each group. This, in turn,
#' allows the estimation of the P and S probabilities, as described by Evett and
#' Buckleton (1990), which used in computing the likelihood ratio (LR) for
#' activity level propositions. The data itself arises from clothing surveys.
#' The zero-inflated zeta distribution has probability mass function
#' \deqn{p(k) = \begin{cases}
#' \pi + \frac{(1-\pi)}{\zeta(s)}&,k=0, \\
#' \frac{(1-\pi)k^{-s}}{\zeta(s)}&,k=1,2,\ldots
#' \end{cases}
#' }{
#' p(0) = pi + (1-pi)/zeta(s),k = 0 and p(k) = (1-pi)k^-s/zeta(s),k>1
#' }
#' where \eqn{\zeta(s)}{zeta(s)} is the Reimann Zeta function.
#'
#' @details The function returns an object of class \code{psFit} which is a
#'   \code{list} contains seven elements:
#' \itemize{
#' \item{\code{psData}}{ -- an object of class \code{psData}--see \code{\link{readData}},}
#' \item{\code{fit}}{ -- the fitted object from \code{\link[stats]{optim}},}
#' \item{\code{pi}}{ - the maximum likelihood estimate of the mixing parameter,}
#' \item{\code{shape}}{ -- the maximum likelihood estimate of the shape parameter,}
#' \item{\code{var.cov}}{ -- the estimated variance-covariance matrix for the parameters,}
#' \item{\code{fitted}}{ -- a named \code{vector} containing the first \code{nterms of
#' the fitted distribution.}}
#' \item{\code{model}}{ -- set to \code{"ziz"} for this model.}
#' }
#'
#' The output can be used in a variety of ways. If the interest is just in the
#' mixing and shape parameter estimates, then the \code{pi} and \code{shape}
#' member of the \code{psFit} object contains this information. It is also
#' displayed along with a number of fitted probabilities by the
#' \code{\link{print.psFit}} method. The fitted object can also be plotted
#' using the plot method \code{\link{plot.psFit}}, and to create a probability
#' function with \code{\link{probfun}}. **NOTE** The value of the shape
#' parameter that is printed (if you print the fitted object) is different
#' from that value that is stored in \code{shape}. The stored value is for the
#' \pkg{VGAM} parameterisation of the zeta distribution which uses
#' \eqn{s^\prime = s - 1}{s' = s - 1}. Therefore the printed value is \eqn{s =
#' s^\prime + 1}{s = s' + 1}. If you intend to use the fitted value with
#' \code{\link[VGAM]{dzeta}}, then you should use the stored value
#' \eqn{s^\prime}{s'}.
#'
#' If \code{start} is not specified, then it is set to (0.5, 1). The reason
#' the starting values are not zero is that small starting values seem to
#' cause instability in the likelihood. If you specify your own starting
#' value, it would be sensible to keep both above 0.5.
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
#' @param start a starting value for the optimiser.
#' @param ... other parameters - not currently used.
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
                     start = c(0.5, 1),
                     ...){
  nvals = 1:nterms
  if(!is(x, "psData")){
    stop("x must be an object of class psData")
  }

  if(start[1] <= 0 || start[1] >= 1){
    stop("The starting value for pi must be in (0, 1)")
  }

  if(start[2] <= 0){
    stop("The Zeta function is undefined for shape = 0. Choose a start value > 0.")
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

  y = rep(obsData, x$data$rn)

  d.one.inflated.zeta = function(x, shape, p, log = FALSE){
    rval = (1 - p) * VGAM::dzeta(x, shape = shape)
    rval[x == 1] = rval[x == 1] + p

    if(log){
      return(log(rval))
    }
    return(rval)
  }

  logLik = function(params){
    p = params[1]
    shape = params[2]

    rval = (1 - p) * VGAM::dzeta(obsData, shape = shape)
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
              lower = c(sqrt(.Machine$double.eps), 0.1),
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

#' @rdname fitZIDist
#' @export
fitZIdist = fitZIDist

#' @rdname fitZIDist
#' @export
fitzidist = fitZIDist

zi.loglik = function(y, theta){
  p = theta[1]
  shape = theta[2]
  rval = (1 - p) * VGAM::dzeta(y$n, shape = shape)
  rval[y$n == 1] = rval[y$n == 1] + p
  sum(y$rn *log(rval))
}


fitZIDistPL = function(x, nterms = 10,
                       start = c(0.5, 1),
                       lambda = 0.1,
                     ...){
  nvals = 1:nterms
  if(!is(x, "psData")){
    stop("x must be an object of class psData")
  }

  if(start[1] <= 0 || start[1] >= 1){
    stop("The starting value for pi must be in (0, 1)")
  }

  if(start[2] <= 0){
    stop("The Zeta function is undefined for shape = 0. Choose a start value > 0.")
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

  y = rep(obsData, x$data$rn)

  d.one.inflated.zeta = function(x, shape, p, log = FALSE){
    rval = (1 - p) * VGAM::dzeta(x, shape = shape)
    rval[x == 1] = rval[x == 1] + p

    if(log){
      return(log(rval))
    }
    return(rval)
  }

  logLik = function(params){
    p = params[1]
    shape = params[2]

    rval = (1 - p) * VGAM::dzeta(obsData, shape = shape)
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
              lower = c(sqrt(.Machine$double.eps), 0.1),
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

