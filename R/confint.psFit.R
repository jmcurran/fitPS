#' S3 confint method for objects of class psFit
#'
#' @param x an object of class \code{psFit}---see fitDist for more details
#' @param level the confidence level required---restricted to [0.75, 1)
#'
#' @return a list with two items: \code{wald} and \code{prof} containing the
#' Wald and profile likelihood confidence intervals respectively for the shape
#' parameter of the fitted Zeta distribution. In general these should be relatively
#' close to each other.
#'
#' @examples
#' data(Psurveys)
#' roux = Psurveys$roux
#' fit = fitDist(roux)
#' confint(fit)
#'
#' @importFrom stats confint
#' @export
confint.psFit = function(x, level = 0.95){
  if(level < 0.75 || level >=1){
    stop("Level must be in the interval [0.75,1)")
  }

  obsData = if(x$psData$type == 'P'){ ## the main difference is that the values need 1 added
    x$psData$data$n + 1
  }else{
    x$psData$data$n
  }

  cstar = qchisq(level, 1)

  profLik = Vectorize(
    function(bd){
      2 * (sum(VGAM::dzeta(rep(obsData, x$psData$data$rn), shape = bd, log = TRUE)) +
         x$fit$value) + cstar
    })

  zstar = qnorm((1 - level) * 0.5, lower.tail = FALSE)
  wald = x$shape + c(-1,1) * zstar * sqrt(x$var.shape)
  xmin = max(1, wald[1] - 0.5 * sqrt(x$var.shape))
  xmax = wald[2] + 0.5 * sqrt(x$var.shape)

  lb = uniroot(profLik, c(xmin, x$shape))$root
  ub = uniroot(profLik, c(x$shape, xmax))$root

  results = list(wald = wald, prof = c(lb, ub))
  results = lapply(results, function(ci){
    names(ci) = paste0(100 * c((1 - level) * 0.5,
                               1 - (1 - level) * 0.5), "%")
    return(ci)
  })
  return(results)
}
