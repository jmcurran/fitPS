#' S3 summary method for an object of class \code{psFit}
#'
#' @param object an object of class \code{psFit}, usually from \code{\link{fitDist}}
#' or \code{\link{fitZIDist}}
#' @param ... other arguments passed to \code{summary}
#'
#' @details
#' Experimental because I am unsure if it is useful. If \code{object} is a zero-inflated zeta fitted object,
#' then the function carries out a likelihood ratio test for the value of pi. Currently not implemented for the logarithmic
#' distribution because we are currently not interested in the logarithmic distribution.
#'
#'
#' @importFrom stats pchisq printCoefmat
#' @return No return value, called for side effects
#' @export
summary.psFit = function(object, ...){
  if(object$model == "zeta"){
    cmat = matrix(c(object$shape, sqrt(object$var.shape)), nrow = 1)
    colnames(cmat) = c("Estimate", "Std.Err")
    rownames(cmat) = "shape"
  }else if(object$model == "ziz"){

    ll.H0 = -fitDist(object$psData)$fit$value
    ll.mle = -object$fit$value
    lrt.Stat = 2 * (ll.mle - ll.H0)

    cmat = rbind(c(object$shape, sqrt(object$var.cov[2,2]), NA, NA),
                 c(object$pi, sqrt(object$var.cov[1,1]), lrt.Stat, pchisq(lrt.Stat, 1, lower.tail = FALSE)))
    colnames(cmat) <- c("Estimate", "Std.Err", "X^2 value", "Pr(>X^2)")
    rownames(cmat) = c("shape", "pi")
  }else{

  }

  printCoefmat(cmat)


}
