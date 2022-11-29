#' S3 plot method for an object of class \code{psFit}
#'
#' @param x an object of class \code{psFit}, usually from \code{\link{fitDist}}
#' @param ylim the limits of the y-axis
#' @param ... other arguments passed to \code{plot}
#'
#' @return No return value, called for side effects
#'
#'
#' @importFrom graphics barplot box legend points
#' @importFrom methods is
#' @export
#'
#' @examples
#' p = readData(system.file("extdata", "p.xlsx", package = "fitPS"))
#' fit = fitDist(p)
#' plot(fit)
plot.psFit = function(x, ylim = c(0, 1), ...){
  if(!is(x, "psFit")){
    stop("x must be an object of class psFit")
  }

  labels = gsub("(P|S)", "", names(x$fitted))

  b = barplot(x$fitted,
              names.arg = labels,
              ylab = "Probability",
              xlab = "n",
              main = if(x$psData$type == 'P'){"P terms"}else{"S terms"},
              ylim = ylim,
              ...)

  phat = x$psData$data$rn / sum(x$psData$data$rn)
  ## got to be a smarter way of doing this
  type = x$psData$type
  obsLabels = paste0(type, x$psData$data$n)
  i = match(obsLabels, paste0(type, labels))
  points(b[i], phat, pch = 'X')
  legend("topright", pch = 'X', legend = "Observed", bty = "n")
  box()
}
