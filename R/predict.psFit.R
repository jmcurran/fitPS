#' S3 predict method for an object of class \code{psFit}
#'
#' @param x an object of class \code{psFit}, usually from \code{\link{fitDist}}
#' @param newdata an optional vector of integers at which to calculate
#' \eqn{\Pr(X = x)}{Pr(X = x)}
#' @param interval either \code{"none"} or \code{"confidence"}, can be
#' abbreviated. If \code{"confidence"} then the bounds of a 100 * \code{level}
#' confidence interval for each predicted probability is provided.
#' @param level the level of a confidence interval. Ignored if
#' \code{interval == "none"}.
#' @param ... other arguments passed to \code{predict}---not used
#'
#' @return either a named vector of fitted probabilities, or a \code{data.frame}
#' with columns \code{predicted}, \code{lower}, and \code{upper} and the row names
#' set to show what terms are being calculatd
#' @examples
#' data(Psurveys)
#' roux = Psurveys$roux
#' fit = fitDist(roux)
#' predict(fit, interval = "conf")
#' @importFrom stats predict
#' @export
predict.psFit = function(x, newdata, interval = c("none", "confidence"), level = 0.95, ...){

  interval = match.arg(interval)
  predicted = NULL

  if(missing(newdata)){
    predicted = x$fitted
    newdata = as.numeric(gsub("^(P|S)([0-9]+)$", "\\2", names(x$fitted)))
  }

  if(any(newdata - floor(newdata) > 0)){
    stop("newdata should only contain integers")
  }

  if(x$psData$type == "S" && any(newdata <= 0)){
    stop("Can only make predictions for size probabilities for values of x >= 1")
  }

  if(is.null(predicted)){
    predicted = VGAM::dzeta(newdata + ifelse(x$psData$type == "P", 1, 0),
                            shape = fit$shape)
  }

  if(interval == "confidence"){
    if(level <= 0.75 || level >= 1){
      stop("Level should be in the interval [0.75, 1)")
    }

    zstar = qnorm((1 - level) * 0.5, lower.tail = FALSE)
    lwr = VGAM::dzeta(newdata + ifelse(x$psData$type == "P", 1, 0),
                      shape = fit$shape - zstar * sqrt(fit$var.shape))
    upr = VGAM::dzeta(newdata + ifelse(x$psData$type == "P", 1, 0),
                      shape = fit$shape + zstar * sqrt(fit$var.shape))

    results = data.frame(predicted = predicted, lower = lwr, upper = upr)
    rownames(results) = paste0(x$psData$type, newdata)

    return(results)
  }else{
    names(predicted) = paste0(x$psData$type, newdata)
    return(predicted)
  }
}

