#' S3 fitted method for an object of class \code{psFit}
#'
#' @param object an object of class \code{psFit}, usually from \code{\link{fitDist}}
#' @param ... other arguments passed to \code{fitted}f---not used
#'
#' @return a named vector of fitted probabilities
#' @importFrom stats fitted
#' @export
fitted.psData = function(object, ...){
  return(object$fitted)
}
