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
#' is the Reimann zeta function. Coulson et al. (2001) did not have an easy way
#' to rapidly compute this quantity, hence their use of linear interpolation.
#'
#' @details The function returns an object of class \code{psFit} which is a \code{list}
#' contains four elements:
#' \itemize{
#' \item{\code{psData}}{ - an object of class \code{psData}--see \code{\link{readData}},}
#' \item{\code{fit}}{ - the fitted object from \code{\link[stats]{nlminb}},}
#' \item{\code{shape}}{ - the maximum likelihood estimate of the shape parameter,}
#' \item{\code{fitted}}{ - a named \code{vector} containing the first \code{nterms of
#' the fitted distribution.}}
#' }.
#' The output can be used in a variety of ways. If the interest is just in the
#' shape parameter estimate, then the \code{shape} member of the \code{psFit} object
#' contains this information. It is also displayed along with a number of
#' fitted probabilities by the \code{\link{print.psFit}} method. The fitted object
#' can also be plotted using the plot method \code{\link{plot.psFit}}, and to
#' create a probability function with \code{\link{probfun}}.
#'
#' @seealso \code{\link{plot.psFit}}, \code{\link{print.psFit}}, \code{\link{probfun}}.
#'
#' @references Coulson, S. A., Buckleton, J. S., Gummer, A. B., and Triggs,
#'   C.M., "Glass on clothing and shoes of members of the general population and
#'   people suspected of breaking crimes", Science & Justice 2001: 41(1): 39--48.
#'
#'   Evett, I. W. and Buckleton, J. S., "The interpretation of glass evidence. A
#'   practical approach", Journal of the Forensic Science Society 1990: 30(4):
#'   215--223.
#'
#' @param x an object of type \code{psData}, usually obtained from
#'   \code{\link{readData}}.
#' @param nterms the number of terms to compute the probability distribution for
#' @param start a starting value for the optimiser
#' @param ... other parameters - not currently used.
#'
#' @importFrom stats nlminb runif
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
                   start = runif(1),
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

  logLik = function(shape){
    -sum(VGAM::dzeta(rep(obsData, x$data$rn), shape = shape, log = TRUE))
  }

  fit = nlminb(start = start,
               objective = logLik,
               lower = 0)

  fitted = VGAM::dzeta(nvals, shape = fit$par)
  names(fitted) = if(x$type == 'P'){
    paste0("P", nvals - 1)
  }else{
    paste0("S", nvals)
  }

  result = list(
    psData = x,
    fit = fit,
    shape = fit$par,
    fitted = fitted
  )

  class(result) = "psFit"

  return(result)
}
