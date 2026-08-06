#' Plot posterior probability summaries
#'
#' Plot posterior means and equal-tailed credible intervals for the fitted P or
#' S probabilities stored in a \code{psPosterior} object.
#'
#' @param x An object of class \code{psPosterior}.
#' @param n \code{NULL}, a single number of leading probabilities to plot, or a
#'   vector of P or S indices to plot.
#' @param showInterval Logical; if \code{TRUE}, draw the stored credible
#'   intervals.
#' @param ylim Optional numeric vector giving the y-axis limits.
#' @param xlab,ylab,main Optional plot labels.
#' @param pch Plotting symbol used for posterior means.
#' @param ... Additional graphical arguments passed to \code{plot}.
#'
#' @return Invisibly returns the plotted posterior probability summary data
#'   frame.
#'
#' @details
#' The points are posterior mean probabilities. The vertical intervals are the
#' equal-tailed credible intervals stored in the posterior object. This method
#' does not change \code{plot.psFit()}, which continues to display fitted model
#' probabilities and observed proportions.
#'
#' @examples
#' data(Psurveys)
#' fit = fitZIDist(
#'   Psurveys$roux,
#'   method = "bayes",
#'   bayesOptions = list(posteriorMethod = "numerical"),
#'   nPiGrid = 31,
#'   nShapeGrid = 31
#' )
#' plot(fit$posterior, n = 5)
#'
#' @importFrom graphics arrows axis plot points
#' @export
plot.psPosterior = function(x,
                             n = NULL,
                             showInterval = TRUE,
                             ylim = NULL,
                             xlab = "Probability term",
                             ylab = "Posterior probability",
                             main = "Posterior probability summaries",
                             pch = 19,
                             ...) {
  if (!inherits(x, "psPosterior")) {
    stop("x must be an object of class psPosterior")
  }

  probabilities = posteriorProbs(x, n = n)
  xPositions = seq_len(nrow(probabilities))

  if (is.null(ylim)) {
    yValues = probabilities$estimate

    if (showInterval) {
      yValues = c(yValues, probabilities$lower, probabilities$upper)
    }

    yValues = yValues[is.finite(yValues)]

    if (length(yValues) == 0L) {
      ylim = c(0, 1)
    } else {
      ylim = range(c(0, yValues))
      if (diff(ylim) == 0) {
        ylim = ylim + c(-0.05, 0.05)
      }
    }
  }

  plot(
    xPositions,
    probabilities$estimate,
    type = "n",
    xaxt = "n",
    xlab = xlab,
    ylab = ylab,
    main = main,
    ylim = ylim,
    ...
  )

  axis(1, at = xPositions, labels = probabilities$term)

  if (showInterval) {
    arrows(
      x0 = xPositions,
      y0 = probabilities$lower,
      x1 = xPositions,
      y1 = probabilities$upper,
      angle = 90,
      code = 3,
      length = 0.05
    )
  }

  points(xPositions, probabilities$estimate, pch = pch)

  invisible(probabilities)
}
