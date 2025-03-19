#' Define a prior density
#'
#' Construct a prior that can be used in the \code{\link{fitdist}} function.
#'
#' @param family One of \code{"loguniform"}, \code{"uniform"} or \code{"custom"}.
#' @param range Optionally the range for which the prior density is
#'              evaluated. It is zero outside of this range. Parameter is
#'              required when \code{family="custom"}.
#' @param logd Optionally (required when \code{family="custom"}.) a function
#'             that evaluates the log density of the prior inside the range.
#'
#' @details The default is a LogUniform[-2,2] prior.
#'
#' @return an object of type \code{psPrior}
#'
#' @export
#'
#' @examples
#' ## With default parameters, the prior will be LogUniform[-2, 2]
#' p1 <- makePrior()
#'
#' # plot the prior density
#' xPlot <- seq(from = 0, to = 10, length = 100)
#' plot(xPlot, exp(p1$logd(xPlot)), type="l")
#'
#' # Alternatively, a Uniform[a, b] prior can be used
#' p2 <- makePrior(family = "uniform", range = c(1,10))
#' plot(xPlot, exp(p2$logd(xPlot)), type="l")
#'
#' # A custom prior needs the log density function to be specified
#' # We define an exponential prior with rate = 1
#' logdexp <- function(x) dexp(x, rate = 1/10, log = TRUE)
#' p3 <- makePrior(family = "custom", range = c(1, 10), logd = logdexp)
#'
#' plot(xPlot, exp(p3$logd(xPlot)), type="l")
#'
#' @seealso readData
makePrior = function(family = c("loguniform", "uniform", "custom"),
                                range, logd){

  family = match.arg(family)

  if (missing(range)){
    if (family == "loguniform"){
      range <- exp(c(-2, 2))
    }else if (family == "uniform"){
      range <- exp(c(-2, 2))
    }else if (family == "custom"){
      stop("range needs to be provided for custom prior")
    }
  }

  if (missing(logd)){
    if (family == "loguniform"){
      logw <- log(range[2] - range[1])

      logd <- function(x) {
        sapply(x, function(x){
          if (!inRange(x, range)) return(-Inf)

          -(logw + log(x))})
        }
    }else if (family == "uniform"){
      range <- c(1, exp(2))

      logw <- log(range[2] - range[1])

      logd <- function(x) ifelse(inRange(x, range), yes = -logw, no = -Inf)

    }else if (family == "custom"){
      stop("range needs to be provided for custom prior")
    }
  }

  prior <- list(family = family, range = range, logd = logd)
  class(prior) <- "psPrior"

  prior
}

inRange = function(x, range){
  (x > range[1]) & (x < range[2])
}
