#' Plot a posterior density for a fitted power-series model
#'
#' Plot the marginal posterior density for a parameter in a Bayesian
#' `psFit` object. The posterior density is estimated from stored MCMC
#' samples when they are available. For numerical integration fits, the stored
#' posterior density function is evaluated on a grid.
#'
#' @param object an object of class \code{psFit}, usually from \code{\link{fitDist}}
#'   or \code{\link{fitZIDist}}.
#' @param parameter character; the posterior parameter to plot. The default is
#'   \code{"shape"}. Zero-inflated Bayesian fits also support \code{"pi"}.
#' @param level numeric; credible level for the interval, if displayed.
#' @param showEstimate logical; if \code{TRUE}, draw a vertical line at the posterior
#'   point estimate stored in the fitted object.
#' @param showInterval logical; if \code{TRUE}, draw vertical lines for the equal-tail
#'   credible interval.
#' @param nGrid integer; number of grid points used when a stored posterior
#'   density function is evaluated directly.
#' @param xlab,ylab,main optional plot labels.
#' @param ... other graphical arguments passed to \code{\link[graphics]{plot}}.
#'
#' @return Invisibly returns a data frame containing the plotted posterior
#'   density grid. The return value has attributes named \code{"estimate"} and
#'   \code{"interval"} when those values are available.
#'
#' @details
#' This function intentionally does not overload \code{\link{plot.psFit}}, which
#' continues to plot fitted probabilities. For MCMC fits, the density is
#' estimated from \code{object$chain}. For numerical integration fits from
#' \code{fitDist(..., method = "integrate")}, the stored posterior density function
#' in \code{object$pdf} is evaluated on an automatically chosen grid.
#'
#' @examples
#' \dontrun{
#' data(Psurveys)
#' roux = Psurveys$roux
#' fit = fitDist(roux, method = "bayes")
#' plotPosterior(fit)
#' }
#'
#' @importFrom graphics abline lines plot polygon
#' @importFrom methods is
#' @importFrom stats approx approxfun density integrate quantile uniroot
#' @export
plotPosterior = function(object,
                         parameter = "shape",
                         level = 0.95,
                         showEstimate = TRUE,
                         showInterval = TRUE,
                         nGrid = 512,
                         xlab = NULL,
                         ylab = "Posterior density",
                         main = NULL,
                         ...){
  validatePosteriorPlotInput(object, parameter, level, nGrid)

  posterior = getPosteriorPlotData(object, parameter, level, nGrid)

  if(is.null(xlab)){
    xlab = parameter
  }

  if(is.null(main)){
    main = paste("Posterior density for", parameter)
  }

  plot(
    posterior$x,
    posterior$density,
    type = "l",
    xlab = xlab,
    ylab = ylab,
    main = main,
    ...
  )

  if(showInterval && !any(is.na(posterior$interval))){
    intervalHeight = approxPosteriorHeight(posterior, posterior$interval)
    polygon(
      c(posterior$interval[1], posterior$interval[1], posterior$interval[2], posterior$interval[2]),
      c(0, intervalHeight[1], intervalHeight[2], 0),
      border = NA,
      density = 20
    )
    lines(posterior$x, posterior$density)
    abline(v = posterior$interval, lty = 2)
  }

  if(showEstimate && is.finite(posterior$estimate)){
    abline(v = posterior$estimate, lty = 3)
  }

  return(invisible(makePosteriorPlotReturn(posterior)))
}

validatePosteriorPlotInput = function(object, parameter, level, nGrid){
  if(!methods::is(object, "psFit")){
    stop("object must be an object of class psFit.", call. = FALSE)
  }

  if(length(object$method) != 1 || !object$method %in% c("bayes", "integrate")){
    stop("plotPosterior is only available for Bayesian psFit objects.", call. = FALSE)
  }

  if(length(parameter) != 1 || !is.character(parameter)){
    stop("parameter must be a single character value.", call. = FALSE)
  }

  if(!parameter %in% c("shape", "pi")){
    stop("parameter must be either \"shape\" or \"pi\".", call. = FALSE)
  }

  if(parameter == "pi" && object$model != "ziz"){
    stop("parameter = \"pi\" is only available for zero-inflated Bayesian fits.", call. = FALSE)
  }

  if(length(level) != 1 || !is.numeric(level) || !is.finite(level) || level <= 0 || level >= 1){
    stop("level must be a single numeric value in (0, 1).", call. = FALSE)
  }

  if(length(nGrid) != 1 || !is.numeric(nGrid) || !is.finite(nGrid) || nGrid < 32){
    stop("nGrid must be a single numeric value greater than or equal to 32.", call. = FALSE)
  }
}

getPosteriorPlotData = function(object, parameter, level, nGrid){
  samples = getPosteriorSamples(object, parameter)

  if(!is.null(samples)){
    return(getSamplePosteriorPlotData(object, parameter, level, samples))
  }

  getDensityPosteriorPlotData(object, parameter, level, nGrid)
}

getPosteriorSamples = function(object, parameter){
  if(is.null(object$chain)){
    return(NULL)
  }

  if(is.data.frame(object$chain) || is.matrix(object$chain)){
    if(!parameter %in% colnames(object$chain)){
      return(NULL)
    }
    samples = object$chain[, parameter]
  }else if(parameter == "shape" && is.numeric(object$chain)){
    samples = object$chain
  }else{
    return(NULL)
  }

  samples = samples[is.finite(samples)]
  if(length(samples) < 2){
    return(NULL)
  }

  samples
}

getSamplePosteriorPlotData = function(object, parameter, level, samples){
  densityEstimate = density(samples)
  alpha = (1 - level) / 2
  interval = as.numeric(quantile(samples, probs = c(alpha, 1 - alpha), names = FALSE))

  list(
    x = densityEstimate$x,
    density = densityEstimate$y,
    estimate = getPosteriorEstimate(object, parameter),
    interval = interval
  )
}

getDensityPosteriorPlotData = function(object, parameter, level, nGrid){
  densityFunction = getPosteriorDensityFunction(object, parameter)

  if(is.null(densityFunction)){
    stop("No posterior samples or stored density function are available for this parameter.", call. = FALSE)
  }

  grid = getPosteriorGrid(object, parameter, nGrid)
  densityValues = pmax(as.numeric(densityFunction(grid)), 0)
  densityValues[!is.finite(densityValues)] = 0

  if(sum(densityValues) <= 0){
    stop("The stored posterior density could not be evaluated on a useful grid.", call. = FALSE)
  }

  posteriorFunctions = makePosteriorGridFunctions(grid, densityValues)
  densityValues = densityValues / posteriorFunctions$area

  alpha = (1 - level) / 2
  interval = posteriorGridQuantile(
    posteriorFunctions$cdf,
    range(grid),
    probs = c(alpha, 1 - alpha)
  )

  list(
    x = grid,
    density = densityValues,
    estimate = getPosteriorEstimate(object, parameter),
    interval = interval
  )
}


getPosteriorDensityFunction = function(object, parameter){
  if(!is.null(object$marginalPdf) && is.list(object$marginalPdf)){
    densityFunction = object$marginalPdf[[parameter]]
    if(is.function(densityFunction)){
      return(densityFunction)
    }
  }

  if(parameter == "shape" && !is.null(object$pdf) && is.function(object$pdf)){
    return(object$pdf)
  }

  NULL
}

getPosteriorGrid = function(object, parameter, nGrid){
  estimate = getPosteriorEstimate(object, parameter)
  spread = getPosteriorSpread(object, parameter)

  if(!is.finite(estimate)){
    estimate = if(parameter == "pi"){
      0.5
    }else{
      2
    }
  }

  if(!is.finite(spread) || spread <= 0){
    spread = max(abs(estimate) / 4, 0.25)
  }

  lower = estimate - 6 * spread
  upper = estimate + 6 * spread

  if(parameter == "shape"){
    lower = max(1 + sqrt(.Machine$double.eps), lower)
    upper = max(upper, lower + 6 * spread)
  }else{
    lower = max(0, lower)
    upper = min(1, upper)
  }

  seq(lower, upper, length.out = as.integer(nGrid))
}

getPosteriorEstimate = function(object, parameter){
  if(!is.null(object[[parameter]]) && is.numeric(object[[parameter]]) && length(object[[parameter]]) == 1){
    return(object[[parameter]])
  }

  if(!is.null(object$fit$par)){
    par = object$fit$par
    if(!is.null(names(par)) && parameter %in% names(par)){
      return(par[[parameter]])
    }

    if(parameter == "shape" && length(par) == 1){
      return(as.numeric(par[1]))
    }
  }

  NA_real_
}

getPosteriorSpread = function(object, parameter){
  if(parameter == "shape" && !is.null(object$var.shape)){
    return(sqrt(as.numeric(object$var.shape)))
  }

  if(!is.null(object$var.cov) && parameter %in% colnames(object$var.cov)){
    return(sqrt(as.numeric(object$var.cov[parameter, parameter])))
  }

  NA_real_
}

makePosteriorGridFunctions = function(x, y){
  densityFun = approxfun(x, y, yleft = 0, yright = 0, rule = 1)
  xRange = range(x)
  area = integrate(densityFun, lower = xRange[1], upper = xRange[2])$value

  if(!is.finite(area) || area <= 0){
    stop("The stored posterior density could not be normalized.", call. = FALSE)
  }

  cdfFun = function(q){
    q = pmin(pmax(q, xRange[1]), xRange[2])
    vapply(
      q,
      function(value){
        integrate(densityFun, lower = xRange[1], upper = value)$value / area
      },
      numeric(1)
    )
  }

  list(
    density = densityFun,
    cdf = cdfFun,
    area = area
  )
}

posteriorGridQuantile = function(cdfFun, xRange, probs){
  vapply(
    probs,
    function(prob){
      if(prob <= 0){
        return(xRange[1])
      }

      if(prob >= 1){
        return(xRange[2])
      }

      uniroot(
        function(x){
          cdfFun(x) - prob
        },
        interval = xRange
      )$root
    },
    numeric(1)
  )
}

approxPosteriorHeight = function(posterior, x){
  approx(posterior$x, posterior$density, xout = x, rule = 2)$y
}

makePosteriorPlotReturn = function(posterior){
  result = data.frame(
    x = posterior$x,
    density = posterior$density
  )
  attr(result, "estimate") = posterior$estimate
  attr(result, "interval") = posterior$interval
  result
}
