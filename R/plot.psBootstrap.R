#' Plot bootstrap probability summaries
#'
#' Plot bootstrap mean probabilities and percentile confidence intervals for
#' the fitted P or S probabilities stored in a \code{psBootstrap} object.
#'
#' @param x An object of class \code{psBootstrap}.
#' @param n \code{NULL}, a single number of leading probabilities to plot, or a
#'   vector of P or S indices to plot.
#' @param showInterval Logical; if \code{TRUE}, draw the stored percentile
#'   confidence intervals.
#' @param ylim Optional numeric vector giving the y-axis limits.
#' @param xlab,ylab,main Optional plot labels.
#' @param pch Plotting symbol used for bootstrap means.
#' @param ... Additional graphical arguments passed to \code{plot}.
#'
#' @return Invisibly returns the plotted bootstrap probability summary data
#'   frame.
#'
#' @examples
#' if (interactive()) {
#' data(Psurveys)
#' fit = fit(
#'   Psurveys$roux,
#'   model = zizModel(),
#'   method = "bootstrap",
#'   nterms = 4,
#'   B = 20,
#'   seed = 123,
#'   silent = TRUE,
#'   parallel = FALSE
#' )
#' plot(fit$bootstrap, n = 4)
#'
#' }
#' @importFrom graphics arrows axis plot points
#' @export
plot.psBootstrap = function(x,
                             n = NULL,
                             showInterval = TRUE,
                             ylim = NULL,
                             xlab = "Probability term",
                             ylab = "Bootstrap probability",
                             main = "Bootstrap probability summaries",
                             pch = 19,
                             ...) {
  if (!inherits(x, "psBootstrap")) {
    stop("x must be an object of class psBootstrap")
  }

  probabilities = bootstrapProbs(x, n = n)
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
