#' Bayesian credible intervals or regions
#'
#' Extract one-dimensional credible intervals or two-dimensional credible
#' regions from a fitted Bayesian uncertainty representation. Parametric
#' Bayesian fits reuse the posterior representation stored by [fit()], while
#' Rubin Bayesian Bootstrap objects reuse their stored weighted-fit parameter
#' replicates.
#'
#' @param psFit A Bayesian `psFit` object or a `psBayesianBootstrap` object.
#' @param level One or more credible probability levels strictly between zero
#'   and one.
#' @param plot Logical; if `TRUE`, draw the corresponding uncertainty display
#'   with [plotUncertainty()] after extracting the interval or region.
#' @param silent Logical; retained for backward compatibility. If `FALSE`, a
#'   short progress message is printed for two-dimensional regions.
#' @param parameters Optional character vector naming one or two parameters.
#'   By default all model parameters are used when there are at most two.
#' @param nGrid Number of grid points used when a one-dimensional posterior
#'   density must be evaluated for interval extraction.
#' @param ... Additional graphical arguments passed to [plotUncertainty()] when
#'   `plot = TRUE`.
#'
#' @aliases credInt
#'
#' @details
#' `credint()` is the Bayesian numerical extractor corresponding to the visual
#' [plotUncertainty()] interface. It does not rerun posterior fitting merely to
#' obtain an interval or region.
#'
#' For one-dimensional numerical posteriors, credible limits are interpolated
#' from the stored cumulative posterior representation. For MCMC posteriors,
#' stored draws are used. Importance-sampling regions preserve their stored
#' importance weights, Laplace regions use the stored Gaussian covariance
#' approximation, and two-dimensional numerical posteriors use the stored
#' posterior grid and quadrature mass. Rubin Bayesian Bootstrap regions are
#' derived from the stored weighted-fit parameter replicates.
#'
#' The returned intervals are equal-tailed in one dimension. Sample-based
#' two-dimensional regions use probability-content KDE contours; numerical
#' posterior grids use highest-density probability-content contours; Laplace
#' regions use the corresponding Gaussian probability ellipse.
#'
#' @return For one parameter, a named numeric vector for one requested level or
#'   a data frame of lower and upper limits for multiple levels. For two
#'   parameters, a named list of credible-region contour coordinates. Each
#'   contour records its probability content and density threshold where one is
#'   available.
#'
#' @examples
#' if (interactive()) {
#'   data(Psurveys)
#'   bayesFit = fit(
#'     Psurveys$roux,
#'     model = zizModel(),
#'     method = "bayes",
#'     bayesOptions = list(posteriorMethod = "numerical")
#'   )
#'   credint(bayesFit, level = c(0.80, 0.95))
#'   plotUncertainty(bayesFit, level = c(0.80, 0.95))
#' }
#'
#' @importFrom stats quantile
#' @export
credint = function(psFit,
                   level = 0.95,
                   plot = FALSE,
                   silent = FALSE,
                   parameters = NULL,
                   nGrid = 401,
                   ...) {
  level = validateUncertaintyLevels(level)
  validateUncertaintyGridSize(nGrid)

  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    stop("plot must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(silent) || length(silent) != 1L || is.na(silent)) {
    stop("silent must be TRUE or FALSE", call. = FALSE)
  }

  uncertainty = if (inherits(psFit, "psBayesianBootstrap")) {
    extractBayesianBootstrapCredibleUncertainty(
      object = psFit,
      level = level,
      parameters = parameters
    )
  } else {
    if (!is(psFit, "psFit")) {
      stop(
        "credint() requires a Bayesian psFit or psBayesianBootstrap object",
        call. = FALSE
      )
    }
    if (!identical(psFit$method, "bayes")) {
      stop(
        "credint() reports Bayesian credible uncertainty; fit the model with method = \"bayes\" first",
        call. = FALSE
      )
    }
    extractParametricBayesianCredibleUncertainty(
      object = psFit,
      level = level,
      parameters = parameters,
      nGrid = nGrid,
      silent = silent
    )
  }

  if (plot) {
    plotUncertainty(
      psFit,
      level = level,
      parameters = uncertainty$parameters,
      nGrid = nGrid,
      ...
    )
  }

  formatCredibleUncertainty(uncertainty)
}

#' Extract credible uncertainty from a parametric Bayesian fit.
#'
#' @keywords internal
#' @noRd
extractParametricBayesianCredibleUncertainty = function(object,
                                                          level,
                                                          parameters,
                                                          nGrid,
                                                          silent) {
  model = modelFromFit(object)
  parameterNames = resolveUncertaintyParameters(
    modelParameterNames(model),
    parameters
  )

  if (length(parameterNames) == 1L) {
    intervals = t(vapply(level, function(probabilityLevel) {
      getPosteriorPlotData(
        object = object,
        parameter = parameterNames[[1L]],
        level = probabilityLevel,
        nGrid = nGrid
      )$interval
    }, numeric(2L)))
    colnames(intervals) = c("lower", "upper")
    rownames(intervals) = as.character(level)

    return(list(
      dimension = 1L,
      parameters = parameterNames,
      level = level,
      intervals = intervals
    ))
  }

  if (!silent) {
    cat("Computing credible regions\n")
  }

  representation = object$posterior$representation
  chain = representation$value$chain
  if (!is.null(chain)) {
    values = as.data.frame(chain)[, parameterNames, drop = FALSE]
    region = makeSampleUncertaintyRegion(values, parameterNames, level)
    return(list(
      dimension = 2L,
      parameters = parameterNames,
      level = level,
      contours = region$contours
    ))
  }

  approximation = representation$value$approximation
  if (!is.null(approximation$samples) &&
      is.data.frame(approximation$samples) &&
      all(c(parameterNames, "weight") %in% names(approximation$samples))) {
    values = approximation$samples[, parameterNames, drop = FALSE]
    region = makeSampleUncertaintyRegion(
      values = values,
      parameters = parameterNames,
      level = level,
      weights = approximation$samples$weight
    )
    return(list(
      dimension = 2L,
      parameters = parameterNames,
      level = level,
      contours = region$contours
    ))
  }

  if (!is.null(approximation$mode) && !is.null(approximation$varCov)) {
    contours = makeGaussianUncertaintyContours(
      centre = approximation$mode[parameterNames],
      covariance = approximation$varCov[
        parameterNames,
        parameterNames,
        drop = FALSE
      ],
      level = level,
      parameters = parameterNames
    )
    return(list(
      dimension = 2L,
      parameters = parameterNames,
      level = level,
      contours = contours
    ))
  }

  grid = representation$value$grid
  if (!is.null(grid) && !is.null(grid$parameters) && !is.null(grid$weights)) {
    region = makeGridUncertaintyRegion(grid, parameterNames, level)
    return(list(
      dimension = 2L,
      parameters = parameterNames,
      level = level,
      contours = region$contours
    ))
  }

  stop(
    "credint() could not extract a supported two-dimensional posterior representation",
    call. = FALSE
  )
}

#' Extract credible uncertainty from Rubin Bayesian Bootstrap replicates.
#'
#' @keywords internal
#' @noRd
extractBayesianBootstrapCredibleUncertainty = function(object,
                                                         level,
                                                         parameters) {
  replicates = object$replicates$parameters
  parameterNames = resolveUncertaintyParameters(names(replicates), parameters)
  values = replicates[, parameterNames, drop = FALSE]
  values = values[complete.cases(values), , drop = FALSE]

  if (nrow(values) < 2L) {
    stop("at least two successful Bayesian Bootstrap replicates are required", call. = FALSE)
  }

  if (length(parameterNames) == 1L) {
    intervals = t(vapply(level, function(probabilityLevel) {
      alpha = (1 - probabilityLevel) / 2
      quantile(
        values[[parameterNames[[1L]]]],
        probs = c(alpha, 1 - alpha),
        names = FALSE,
        type = 7
      )
    }, numeric(2L)))
    colnames(intervals) = c("lower", "upper")
    rownames(intervals) = as.character(level)

    return(list(
      dimension = 1L,
      parameters = parameterNames,
      level = level,
      intervals = intervals
    ))
  }

  region = makeSampleUncertaintyRegion(values, parameterNames, level)
  list(
    dimension = 2L,
    parameters = parameterNames,
    level = level,
    contours = region$contours
  )
}

#' Format a common credible-uncertainty representation for public return.
#'
#' @keywords internal
#' @noRd
formatCredibleUncertainty = function(uncertainty) {
  if (identical(uncertainty$dimension, 1L)) {
    intervals = uncertainty$intervals
    if (nrow(intervals) == 1L) {
      result = as.numeric(intervals[1L, ])
      names(result) = c("lower", "upper")
      return(result)
    }
    return(as.data.frame(intervals))
  }

  parameters = uncertainty$parameters
  levelNames = paste0(format(100 * uncertainty$level, trim = TRUE), "%")
  regions = lapply(uncertainty$level, function(probabilityLevel) {
    matchingContours = Filter(
      function(contour) isTRUE(all.equal(contour$level, probabilityLevel)),
      uncertainty$contours
    )
    lapply(matchingContours, function(contour) {
      result = list(
        level = contour$densityLevel,
        probability = contour$level
      )
      result[[parameters[[1L]]]] = contour$x
      result[[parameters[[2L]]]] = contour$y
      result
    })
  })
  names(regions) = levelNames

  if (length(regions) > 1L && all(vapply(regions, length, integer(1L)) == 1L)) {
    regions = lapply(regions, `[[`, 1L)
  }
  regions
}

#' @describeIn credint Bayesian credible intervals or regions
#' @export
credInt = credint
