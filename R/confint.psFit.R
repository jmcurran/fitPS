#' S3 confint method for objects of class psFit
#'
#' @param object an object of class \code{psFit}---see fitDist for more details
#' @param parm added for compatibility. Should be left empty as it is ignored.
#' @param level the confidence level required---restricted to [0.75, 1)
#' @param ... in theory other parameters to be passed to \code{confint}, but in
#'   reality ignored.
#'
#' @details NOTE: this method only works for the Zeta model. If you need a
#'   confidence region for the ZIZ model, then use bootstrapping via the
#'   \code{\link{bootCI}} function.
#'
#' @return a list with two items: \code{wald} and \code{prof} containing the
#'   Wald and profile likelihood confidence intervals respectively for the shape
#'   parameter of the fitted Zeta distribution. In general these should be
#'   relatively close to each other. **NOTE** These values are for the
#'   \pkg{VGAM} parameterisation of the Zeta distribution which uses
#'   \eqn{s^\prime = s - 1}{s' = s - 1}. This means they can be used without
#'   alteration in \code{\link[VGAM]{dzeta}}.
#'
#'
#' @examples
#' data(Psurveys)
#' roux = Psurveys$roux
#' fit = fitDist(roux)
#' confint(fit)
#'
#' @importFrom stats confint qchisq uniroot
#' @export
confint.psFit = function(object, parm, level = 0.95, ...){
  if(object$zeroInflated){
    msg = paste0("This method does not currently work for the ZIZ model. However",
                  "the bootCI function will give you a confidence region",
                 collapse = "\n")
    stop(msg)
  }

  if(level < 0.75 || level >=1){
    stop("Level must be in the interval [0.75,1)")
  }

  obsData = if(object$psData$type == 'P'){ ## the main difference is that the values need 1 added
    object$psData$data$n + 1
  }else{
    object$psData$data$n
  }

  cstar = qchisq(level, 1)

  profLik = Vectorize(
    function(bd){
      2 * (sum(VGAM::dzeta(rep(obsData, object$psData$data$rn), shape = bd, log = TRUE)) +
         object$fit$value) + cstar
    })

  zstar = qnorm((1 - level) * 0.5, lower.tail = FALSE)
  wald = object$shape + c(-1,1) * zstar * sqrt(object$var.shape)
  xmin = max(0, wald[1] - 0.5 * sqrt(object$var.shape))
  xmax = wald[2] + 0.5 * sqrt(object$var.shape)

  lb = uniroot(profLik, c(xmin, object$shape), extendInt = "yes")$root
  ub = uniroot(profLik, c(object$shape, xmax), extendInt = "yes")$root

  results = list(wald = wald, prof = c(lb, ub))
  results = lapply(results, function(ci){
    names(ci) = paste0(100 * c((1 - level) * 0.5,
                               1 - (1 - level) * 0.5), "%")
    return(ci)
  })
  return(results)
}
