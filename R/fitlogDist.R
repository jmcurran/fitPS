#' Fit a Logarithmic Distribution to Forensic Data
#'
#' This function uses maximum likelihood estimation (MLE) to estimate the shape
#' parameter of a logarithmic distribution from a set of observed counts
#' for either the number of groups/sources of forensically interesting material
#' (mostly glass or paint) recovered from clothing, or the number of
#' fragments/particles in each group. This, in turn, allows the estimation of
#' the P and S probabilities, as described by Evett and Buckleton (1990), which
#' used in computing the likelihood ratio (LR) for activity level propositions.
#' The data itself arises from clothing surveys. The logarithmic
#' distribution has probability mass function
#' \deqn{p(k) = \frac{\pi^k}{k\log_e(1-\pi)},0<\pi<1.}
#' {pi^k/(k * ln(1-pi)), 0<pi<1.}
#'
#' @details The function returns an object of class \code{psFit} which is a
#'   \code{list} contains seven elements:
#' \itemize{
#' \item{\code{psData}}{ -- an object of class \code{psData}--see \code{\link{readData}},}
#' \item{\code{fit}}{ -- the fitted object from \code{\link[stats]{optim}},}
#' \item{\code{pi}}{ - the maximum likelihood estimate of the shape parameter,}
#' \item{\code{var}}{ -- the estimated variance for the shape parameter,}
#' \item{\code{fitted}}{ -- a named \code{vector} containing the first \code{nterms of
#' the fitted distribution.}}
#' }
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
#' @keywords internal
#' @export
#'
#' @examples
#' data(Psurveys)
#' roux = Psurveys$roux
#' fit = fitlogDist(roux)
#' fit
fitlogDist = function(x, nterms = 10,
                     start = 0.5,
                     ...){
  nvals = 1:nterms
  if(!is(x, "psData")){
    stop("x must be an object of class psData")
  }

  if(start <= 0 || start >= 1){
    stop("The starting value for pi must be in (0, 1)")
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

  logLik = function(params){
    p = params

    rval = sum(x$data$rn * VGAM::dlog(obsData, shape = p, log = TRUE))

    if(is.infinite(rval) || is.nan(rval)){
      stop(sprintf("Infinite log-likelihod: pi = %6.4E\n", p))
    }
    return(-rval)
  }

  # fit = nlminb(start = start,
  #              objective = logLik,
  #              lower = 1)

  fit = optim(par = start,
              fn = logLik,
              method = "L-BFGS-B",
              lower = c(sqrt(.Machine$double.eps)),
              upper  = c(1 - .Machine$double.eps),
              hessian = TRUE)

  fitted = VGAM::dlog(nvals, shape = fit$par)
  names(fitted) = if(x$type == 'P'){
    paste0("P", nvals - 1)
  }else{
    paste0("S", nvals)
  }

  result = list(
    psData = x,
    fit = fit,
    pi = fit$par,
    var = 1 / fit$hessian[1,1],
    fitted = fitted,
    model = "log"
  )


  class(result) = "psFit"

  return(result)
}

#' @rdname fitlogDist
#' @export
fitLogdist = fitlogDist

#' @rdname fitlogDist
#' @export
fitlogdist = fitlogDist

