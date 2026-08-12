#' Plug-in probability functions
#'
#' Creates a probability function that evaluates P or S terms at the fitted
#' parameter values stored in a \code{psFit} object.
#'
#' @param psFitobj an object of class \code{psFit}--see \code{\link{fitDist}}
#' and \code{\link{fitZIDist}}.
#'
#' @details
#' This is deliberately a plug-in probability function. For Bayesian fits,
#' use \code{\link{posteriorProbs}} for posterior mean probabilities. For
#' maximum-likelihood fits with an attached bootstrap distribution, use
#' \code{\link{bootstrapProbs}} for bootstrap mean probabilities and percentile
#' confidence intervals.
#'
#' @return a function that can be used to calculate any plug-in P or S term.
#' @export
#'
#' @examples
#' p = readData(system.file("extdata", "p.xlsx", package = "fitPS"))
#' fit = fitDist(p)
#' P = probfun(fit)
#' P(0:5)
probfun = function(psFitobj) {
  pf = function(x) {
    if (psFitobj$model == "zeta") {
      if (psFitobj$psData$type == "P") {
        p = dzetaStandard(x + 1, shape = psFitobj$shape)
        names(p) = paste0("P", x)
        return(p)
      }
      p = dzetaStandard(x, shape = psFitobj$shape)
      names(p) = paste0("S", x)
      return(p)
    }

    if (psFitobj$model == "ziz") {
      if (psFitobj$psData$type == "P") {
        p = (1 - psFitobj$pi) *
          dzetaStandard(x + 1, shape = psFitobj$shape)
        p[x == 0] = p[x == 0] + psFitobj$pi
        names(p) = paste0("P", x)
        return(p)
      }
      p = (1 - psFitobj$pi) * dzetaStandard(x, shape = psFitobj$shape)
      p[x == 1] = p[x == 1] + psFitobj$pi
      names(p) = paste0("S", x)
      return(p)
    }

    stop("This function is not currently implemented for the logarithmic distribution.")
  }
  pf
}
