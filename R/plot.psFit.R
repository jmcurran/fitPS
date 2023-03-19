#' S3 plot method for an object of class \code{psFit}
#'
#' @param x an object of class \code{psFit}, usually from \code{\link{fitDist}}
#' @param ylim the limits of the y-axis
#' @param conf if \code{TRUE}, and the model is the the Zeta model (as opposed
#' to the Zero-Inflated Zeta (ZIZ), then confidence intervals (based on the standard
#'   error of the shape parameter) are drawn on the plot. If the ZIZ model has
#'   been used, then this is ignored.
#' @param conf.level the confidence level for the confidence intervals. Must be
#'   between 0.75 and 0.99.
#' @param ci.type Specifies the type of confidence interval. If \code{conf ==
#'   TRUE}, then then \code{ci.type} can be either \code{"wald"} \code{"prof"} (or
#'   an abbreviation), depending on whether the Wald interval or the profile
#'   likelihood interval should be used. Note that these are intervals on the shape
#'   parameter and not the density heights. Therefore the intervals around the
#'   probabilities should not really be thought of as confidence intervals but
#'   rather something more similar to a "sensitivity" interval.
#' @param log.scale if \code{TRUE} the \eqn{y}{y}-axis is changed to a
#'   logarithmic (base 10) axis.
#' @param ... other arguments passed to \code{plot}.
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
#'
#' ## An example with Wald generated intervals
#' plot(fit, conf = TRUE)
#'
## An example with profile likelihood generated intervals
#' plot(fit, conf = TRUE, ci.type = "p")
plot.psFit = function(x, ylim = c(0, 1), conf = FALSE, conf.level = 0.95,
                      ci.type = c("wald", "prof"),
                      log.scale = FALSE, ...){
  if(!is(x, "psFit")){
    stop("x must be an object of class psFit")
  }

  labels = gsub("(P|S)", "", names(x$fitted))

  if(!log.scale){
    b = barplot(x$fitted,
                names.arg = labels,
                ylab = "Probability",
                xlab = "n",
                main = if(x$psData$type == 'P'){"P terms"}else{"S terms"},
                ...)

    phat = x$psData$data$rn / sum(x$psData$data$rn)
    ## got to be a smarter way of doing this
    type = x$psData$type
    obsLabels = paste0(type, x$psData$data$n)
    i = match(obsLabels, paste0(type, labels))
    points(b[i], phat, pch = 'X')
    legend("topright", pch = 'X', legend = "Observed", bty = "n")
    box()

    if(conf && !x$zeroInflated){
      ci.type = match.arg(ci.type)

      if(conf.level < 0.75 || conf.level > 0.99){
        stop("The confidence level must be a sensible value, i.e in [0.75, 0.99].")
      }

      ci = confint(x, level = conf.level)[[ci.type]]

      lbShape = ci[1]
      ubShape = ci[2]
      nvals = 1:length(b)
      lb = VGAM::dzeta(nvals, lbShape)
      ub = VGAM::dzeta(nvals, ubShape)
      Hmisc::errbar(b, x$fitted, ub, lb, add = TRUE)
    }
  }else{
    b = barplot(log10(x$fitted),
                names.arg = labels,
                ylab = expression(log[10](Probability)),
                xlab = "n",
                main = if(x$psData$type == 'P'){"P terms"}else{"S terms"},
                ylim = ylim,
                ...)

    phat = log10(x$psData$data$rn / sum(x$psData$data$rn))
    ## got to be a smarter way of doing this
    type = x$psData$type
    obsLabels = paste0(type, x$psData$data$n)
    i = match(obsLabels, paste0(type, labels))
    points(b[i], phat, pch = 'X')
    legend("topright", pch = 'X', legend = "Observed", bty = "n")
    box()

    if(conf && !x$zeroInflated){
      ci.type = match.arg(ci.type)

      if(conf.level < 0.75 || conf.level > 0.99){
        stop("The confidence level must be a sensible value, i.e in [0.75, 0.99].")
      }

      ci = confint(x, level = conf.level)[[ci.type]]

      lbShape = ci[1]
      ubShape = ci[2]
      nvals = 1:length(b)
      lb = log10(VGAM::dzeta(nvals, lbShape))
      ub = log10(VGAM::dzeta(nvals, ubShape))
      Hmisc::errbar(b, log10(x$fitted), ub, lb, add = TRUE)
    }
  }
}
