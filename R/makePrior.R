#' Define a prior density
#'
#' Construct a prior that can be used in the \code{\link{fitdist}} function.
#'
#' @param family One of \code{"loguniform"}, \code{"uniform"} or \code{"custom"}.
#' @param range Optionally the range for which the prior density is
#'              evaluated. It is zero outside of this range. For zeta models
#'              this range is on the standard \code{shape} scale, where
#'              \code{shape > 1}.
#' @param logd Optionally (required when \code{family="custom"}.) a function
#'             that evaluates the log density of the prior inside the range.
#'
#' @details The default is a LogUniform[-2, 2] prior on \code{shape - 1}.
#' Equivalently, the prior range on the fitPS standard zeta \code{shape}
#' parameter is \code{1 + exp(c(-2, 2))}.
#'
#' @return an object of type \code{psPrior}
#'
#' @export
#'
#' @examples
#' ## With default parameters, the prior will be LogUniform[-2, 2]
#' ## on shape - 1, so the prior support is above shape = 1.
#' p1 <- makePrior()
#'
#' # plot the prior density
#' xPlot <- seq(from = 1.01, to = 10, length = 100)
#' plot(xPlot, exp(p1$logd(xPlot)), type = "l")
#'
#' # Alternatively, a Uniform[a, b] prior can be used on standard shape
#' p2 <- makePrior(family = "uniform", range = c(1.01, 10))
#' plot(xPlot, exp(p2$logd(xPlot)), type = "l")
#'
#' # A custom prior needs the log density function to be specified
#' # We define an exponential prior with rate = 1 on shape - 1
#' logdexp <- function(x) dexp(x - 1, rate = 1 / 10, log = TRUE)
#' p3 <- makePrior(family = "custom", range = c(1.01, 10), logd = logdexp)
#'
#' plot(xPlot, exp(p3$logd(xPlot)), type = "l")
#'
#' @seealso readData
makePrior = function(family = c("loguniform", "uniform", "custom"),
                     range,
                     logd) {
  family = match.arg(family)

  if (missing(range)) {
    range = 1 + exp(c(-2, 2))
  }

  validatePriorRange(range)

  if (missing(logd)) {
    if (family == "loguniform") {
      logWidth = log(range[2] - range[1])

      logd = function(x) {
        sapply(x, function(value) {
          if (!inRange(value, range)) {
            return(-Inf)
          }

          shiftedValue = value - 1
          if (shiftedValue <= 0) {
            return(-Inf)
          }

          -(logWidth + log(shiftedValue))
        })
      }
    } else if (family == "uniform") {
      logWidth = log(range[2] - range[1])

      logd = function(x) {
        ifelse(inRange(x, range), yes = -logWidth, no = -Inf)
      }
    } else if (family == "custom") {
      stop("logd needs to be provided for custom prior")
    }
  }

  prior = list(family = family, range = range, logd = logd)
  class(prior) = "psPrior"

  prior
}

validatePriorRange = function(range) {
  if (!is.numeric(range) || length(range) != 2L || any(!is.finite(range))) {
    stop("range must be a numeric vector of length two")
  }

  if (range[1] <= 1) {
    stop("prior range must be on the standard shape scale with lower bound greater than 1")
  }

  if (range[2] <= range[1]) {
    stop("prior range upper bound must be greater than lower bound")
  }
}

inRange = function(x, range) {
  (x > range[1]) & (x < range[2])
}
