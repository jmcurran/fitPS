#' Fit a Zeta Distribution to Forensic Data
#'
#' This function uses maximum likelihood estimation (MLE) to estimate the shape
#' parameter of  a zeta distribution from a set of observed counts for either
#' the number of groups/sources of forensically interesting material (mostly
#' glass or paint) recovered from clothing, or the number of fragments/particles
#' in each group. This, in turn, allows the estimation of the P and S
#' probabilities, as described by Evett and Buckleton (1990), which used in
#' computing the likelihood ratio (LR) for activity level propositions. The data
#' itself arises from clothing surveys. The general method is described in
#' Coulson et al. (2001), although poor typesetting, and a lack of definition of
#' terms makes it hard to see. This package improves on the estimation in that
#' linear interpolation is not required, and standard numerical optimisation is
#' used instead. The zeta distribution has probability mass function \deqn{p(k)
#' = \frac{k^{-s}}{\zeta(s)}}{p(k) = k^-s/zeta(s)} where \eqn{\zeta(s)}{zeta(s)}
#' is the Reimann Zeta function. Coulson et al. (2001) did not have an easy way
#' to rapidly compute this quantity, hence their use of linear interpolation.
#'
#' @aliases fitdist
#'
#' @details The function returns an object of class \code{psFit} which is a
#'   \code{list} contains four elements:
#' \describe{
#' \item{\code{psData}}{ -- an object of class \code{psData}--see \code{\link{readData}},}
#' \item{\code{fit}}{ -- the fitted object from \code{\link[stats]{optim}},}
#' \item{\code{shape}}{ -- the maximum likelihood estimate of the shape parameter,}
#' \item{\code{var.shape}}{ - the maximum likelihood estimate of the shape parameter,}
#' \item{\code{fitted}}{ - a named \code{vector} containing the first \code{nterms of
#' the fitted distribution.}}
#' \item{\code{model}}{ - set to \code{"zeta"} for this model}
#' }.
#' The output can be used in a variety of ways. If the interest is just in the
#' shape parameter estimate, then the \code{shape} member of the \code{psFit}
#' object contains this information. It is also displayed along with a number
#' of fitted probabilities by the \code{\link{print.psFit}} method. The fitted
#' object can also be plotted using the plot method \code{\link{plot.psFit}},
#' and to create a probability function with \code{\link{probfun}}. **NOTE**
#' The value of the shape parameter that is printed (if you print the fitted
#' object) is different from that value that is stored in \code{shape}. The
#' stored value is for the \pkg{VGAM} parameterisation of the Zeta
#' distribution which uses \eqn{s^\prime = s - 1}{s' = s - 1}. Therefore the
#' printed value is \eqn{s = s^\prime + 1}{s = s' + 1}. If you intend to use
#' the fitted value with \code{\link[VGAM]{dzeta}}, then you should use the
#' stored value \eqn{s^\prime}{s'}.
#'
#' If \code{start} is not specified, then it is chosen randomly from (0.5, 1).
#' The reason the lower value is not zero is that small starting values seem
#' to cause instability in the likelihood. If you specify your own starting
#' value, it would be sensible to keep it above 0.5.
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
#' p = readData(system.file("extdata", "p.xlsx", package = "fitPS"))
#' fit = fitDist(p)
#' fit
fitDist = function(x, nterms = 10,
                   start = 1,
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

  if(start <= 0){
    msg = paste0("The Zeta function is undefined for shape = 1.",
                 "This software uses s = shape - 1",
                 "Choose a start value > 0.", collapse = "\n")
    stop(msg)
  }

  obsData = if(x$type == 'P'){ ## the main difference is that the values need 1 added
                x$data$n + 1
            }else{
              x$data$n
            }

  logLik = function(shape){
    #-sum(VGAM::dzeta(rep(obsData, x$data$rn), shape = shape, log = TRUE))
    -sum(x$data$rn * VGAM::dzeta(obsData, shape = shape, log = TRUE))
  }

  # fit = nlminb(start = start,
  #              objective = logLik,
  #              lower = 1)

  fit = optim(par = start,
                fn = logLik,
                method = "L-BFGS-B",
                lower = 0,
                hessian = TRUE)

  shape = fit$par ## NOTE VGAM's dzeta is parameterised in terms of
                  ## s = alpha - 1. This has consequences in the formula below
  N = sum(x$data$rn)

  vshape = function(shape){
    z = VGAM::zeta(shape + 1)
    zprime = VGAM::zeta(shape + 1, 1)
    zprimeprime = VGAM::zeta(shape + 1, 2)
    numer = z^2
    denom = N *(z * zprimeprime - zprime^2)

    return(numer / denom)
  }

  var.shape = vshape(shape)

  fitted = VGAM::dzeta(nvals, shape = shape)
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
    model = "zeta"
  )

  class(result) = "psFit"

  return(result)
}

#' @describeIn fitDist Fit a Zeta Distribution to Forensic Data
#' export
fitdist = fitDist
