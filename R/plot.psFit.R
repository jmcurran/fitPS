#' S3 plot method for an object of class \code{psFit}
#'
#' @param x an object of class \code{psFit}, usually from \code{\link{fitDist}}
#' @param ylim the limits of the y-axis
#' @param conf if \code{TRUE}, then confidence intervals (based on the standard
#' error of the shape parameter) are drawn on the plot
#' @param conf.level the confidence level for the confidence intervals. Must be
#' between 0.75 and 0.99.
#' @param ... other arguments passed to \code{plot}
#'
#' @return No return value, called for side effects
#'
#'
#' @importFrom graphics barplot box legend points
#' @importFrom Hmisc errbar
#' @importFrom methods is
#' @importFrom stats qnorm
#' @export
#'
#' @examples
#' p = readData(system.file("extdata", "p.xlsx", package = "fitPS"))
#' fit = fitDist(p)
#' plot(fit)
plot.psFit = function(x, ylim = c(0, 1), conf = FALSE, conf.level = 0.95, ...){
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

  if(conf){
    if(conf.level < 0.75 || conf.level > 0.99){
      stop("The confidence level must be a sensible value, i.e in [0.75, 0.99].")
    }

    z = -qnorm((1 - conf.level))
    lbShape = max(x$shape - z * sqrt(x$var.shape), 0)
    ubShape = x$shape + z * sqrt(x$var.shape)
    nvals = 1:length(b)
    lb = VGAM::dzeta(nvals, lbShape)
    ub = VGAM::dzeta(nvals, ubShape)
    Hmisc::errbar(b, x$fitted, ub, lb, add = TRUE)
  }
}
