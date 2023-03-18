#' S3 predict method for an object of class \code{psFit}
#'
#' @param object an object of class \code{psFit}, usually from \code{\link{fitDist}}
#' @param newdata an optional vector of integers at which to calculate
#'   \eqn{\Pr(X = x)}{Pr(X = x)}
#' @param interval either \code{"none"}, \code{"prof"}, or \code{"wald"} and can
#'   be abbreviated. If \code{"prof"} or \code{"wald"}  then an interval, based
#'   on the bounds of a 100 * \code{level} confidence interval for the shape
#'   parameter, is given for each predicted probability. The interval is provided
#'   based on either a Profile Likelihood, or a Wald, confidence interval
#'   for the shape, and therefore cannot really be regarded as a confidence
#'   interval for the probabilities. The intervals might be more sensibly regarded
#'   as a measure of how sensitive the probabilities are to the choice of shape
#'   parameter.
#' @param level the level of a confidence interval. Ignored if \code{interval ==
#'   "none"}.
#' @param ... other arguments passed to \code{predict}---not used
#'
#' @return either a named vector of fitted probabilities, or a \code{data.frame}
#'   with columns \code{predicted}, \code{lower}, and \code{upper} and the row
#'   names set to show what terms are being calculatd
#' @examples
#' data(Psurveys)
#' roux = Psurveys$roux
#' fit = fitDist(roux)
#' predict(fit, interval = "prof")
#' @importFrom stats predict
#' @export
predict.psFit = function(object, newdata, interval = c("none", "prof", "wald"),
                         level = 0.95, ...){

  interval = match.arg(interval)
  predicted = NULL

  if(missing(newdata)){
    predicted = object$fitted
    newdata = as.numeric(gsub("^(P|S)([0-9]+)$", "\\2", names(object$fitted)))
  }

  if(any(newdata - floor(newdata) > 0)){
    stop("newdata should only contain integers")
  }

  if(object$psData$type == "S" && any(newdata <= 0)){
    stop("Can only make predictions for size probabilities for values of n >= 1")
  }

  if(is.null(predicted)){
    predicted = VGAM::dzeta(newdata + ifelse(object$psData$type == "P", 1, 0),
                            shape = object$shape)
  }

  if(interval %in% c("prof", "wald")){
    if(level <= 0.75 || level >= 1){
      stop("Level should be in the interval [0.75, 1)")
    }

    zstar = qnorm((1 - level) * 0.5, lower.tail = FALSE)
    lwr = VGAM::dzeta(newdata + ifelse(object$psData$type == "P", 1, 0),
                      shape = object$shape - zstar * sqrt(object$var.shape))
    upr = VGAM::dzeta(newdata + ifelse(object$psData$type == "P", 1, 0),
                      shape = object$shape + zstar * sqrt(object$var.shape))

    results = data.frame(predicted = predicted, lower = lwr, upper = upr)
    rownames(results) = paste0(object$psData$type, newdata)

    return(results)
  }else{
    names(predicted) = paste0(object$psData$type, newdata)
    return(predicted)
  }
}

