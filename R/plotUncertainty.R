#' Plot parameter uncertainty for a fitPS inferential result
#'
#' Draw a common one- or two-parameter uncertainty display while preserving the
#' statistical interpretation of the fitted method. Maximum-likelihood fits use
#' profile-likelihood confidence intervals or regions, ordinary bootstrap fits use
#' the stored bootstrap parameter distribution, parametric Bayesian fits use the
#' stored posterior representation, and Rubin Bayesian Bootstrap fits use the
#' stored weighted-fit parameter distribution.
#'
#' @param object A `psFit`, `psBootstrap`, or `psBayesianBootstrap` object.
#' @param level One or more probability levels strictly between zero and one.
#' @param parameters Optional character vector naming one or two parameters to
#'   display. By default all available parameters are used when there are at
#'   most two.
#' @param showPoints Logical or `NULL`; for sample-based two-dimensional displays,
#'   show the stored parameter replicates behind the contours. When `NULL`,
#'   ordinary and Rubin Bayesian Bootstrap displays show their realizations by
#'   default, while parametric Bayesian displays do not.
#' @param nGrid Number of points used for one-dimensional likelihood/density
#'   displays.
#' @param xlab,ylab,main Optional plot labels.
#' @param ... Additional graphical arguments passed to the initial plot call.
#'
#' @return Invisibly returns a list describing the plotted uncertainty
#'   representation, including its inferential interpretation and interval or
#'   contour coordinates. Sample-based two-dimensional results also include
#'   empirical containment of the stored realizations for each contour.
#'
#' @details
#' The visual grammar is deliberately consistent across inferential methods,
#' but the regions are not interchangeable. MLE regions are profile-likelihood
#' confidence regions. Ordinary bootstrap regions are smoothed confidence
#' regions based on the bootstrap estimator distribution. Parametric Bayesian
#' regions are posterior credible regions. Rubin Bayesian Bootstrap regions
#' describe the distribution of weighted maximum-likelihood fits induced by
#' Bayesian Bootstrap weights.
#'
#' For two-dimensional numerical Bayesian fits, fitPS reuses the stored
#' posterior summary grid. Probability-content contours are obtained by ordering the stored grid points
#' by posterior density and cumulatively summing their stored quadrature mass;
#' plotting therefore does not rerun adaptive cubature. For MCMC and bootstrap
#' representations, contours are obtained from the stored parameter draws using
#' kernel density estimation. Ordinary and Rubin Bayesian Bootstrap contours use
#' an unconstrained bivariate KDE so that boundary transformations do not distort
#' the observed joint replicate geometry. Importance-sampling posteriors retain
#' their sample weights for weighted KDE contours, while Laplace fits use Gaussian contours
#' implied by the stored local covariance approximation.
#'
#' @examples
#' if (interactive()) {
#'   data(Psurveys)
#'   roux = Psurveys$roux
#'
#'   mleFit = fit(roux, model = zizModel())
#'   plotUncertainty(mleFit, level = c(0.80, 0.95))
#'
#'   bayesFit = fit(
#'     roux,
#'     model = zizModel(),
#'     method = "bayes",
#'     bayesOptions = list(posteriorMethod = "numerical")
#'   )
#'   plotUncertainty(bayesFit, level = c(0.80, 0.95))
#' }
#'
#' @importFrom grDevices contourLines
#' @importFrom graphics abline lines plot points polygon
#' @importFrom ks contourLevels Hscv kde
#' @importFrom stats density qchisq quantile
#' @export
plotUncertainty = function(object, ...) {
  UseMethod("plotUncertainty")
}

#' @rdname plotUncertainty
#' @export
plotUncertainty.psFit = function(object,
                                  level = c(0.80, 0.95),
                                  parameters = NULL,
                                  showPoints = NULL,
                                  nGrid = 401,
                                  xlab = NULL,
                                  ylab = NULL,
                                  main = NULL,
                                  ...) {
  validateUncertaintyLevels(level)
  validateUncertaintyGridSize(nGrid)

  model = modelFromFit(object)
  parameterNames = resolveUncertaintyParameters(
    modelParameterNames(model),
    parameters
  )

  if (identical(object$method, "bootstrap")) {
    showPoints = resolveUncertaintyShowPoints(showPoints, default = TRUE)
    if (!inherits(object$bootstrap, "psBootstrap")) {
      stop("bootstrap psFit object does not contain a psBootstrap distribution", call. = FALSE)
    }
    return(plotUncertaintyReplicates(
      replicates = object$bootstrap$replicates,
      level = level,
      parameters = parameterNames,
      interpretation = "Frequentist bootstrap confidence",
      showPoints = showPoints,
      xlab = xlab,
      ylab = ylab,
      main = main,
      ...
    ))
  }

  if (identical(object$method, "bayes")) {
    showPoints = resolveUncertaintyShowPoints(showPoints, default = FALSE)
    return(plotUncertaintyBayes(
      object = object,
      level = level,
      parameters = parameterNames,
      showPoints = showPoints,
      nGrid = nGrid,
      xlab = xlab,
      ylab = ylab,
      main = main,
      ...
    ))
  }

  if (!identical(object$method, "mle")) {
    stop("plotUncertainty() requires an MLE, bootstrap, or Bayesian psFit object", call. = FALSE)
  }

  plotUncertaintyMle(
    object = object,
    model = model,
    level = level,
    parameters = parameterNames,
    nGrid = nGrid,
    xlab = xlab,
    ylab = ylab,
    main = main,
    ...
  )
}

#' @rdname plotUncertainty
#' @export
plotUncertainty.psBootstrap = function(object,
                                        level = c(0.80, 0.95),
                                        parameters = NULL,
                                        showPoints = TRUE,
                                        nGrid = 401,
                                        xlab = NULL,
                                        ylab = NULL,
                                        main = NULL,
                                        ...) {
  validateUncertaintyLevels(level)
  validateUncertaintyGridSize(nGrid)
  parameterNames = resolveUncertaintyParameters(names(object$replicates), parameters)
  showPoints = resolveUncertaintyShowPoints(showPoints, default = TRUE)

  plotUncertaintyReplicates(
    replicates = object$replicates,
    level = level,
    parameters = parameterNames,
    interpretation = "Frequentist bootstrap confidence",
    showPoints = showPoints,
    xlab = xlab,
    ylab = ylab,
    main = main,
    ...
  )
}

#' @rdname plotUncertainty
#' @export
plotUncertainty.psBayesianBootstrap = function(object,
                                                level = c(0.80, 0.95),
                                                parameters = NULL,
                                                showPoints = TRUE,
                                                nGrid = 401,
                                                xlab = NULL,
                                                ylab = NULL,
                                                main = NULL,
                                                ...) {
  validateUncertaintyLevels(level)
  validateUncertaintyGridSize(nGrid)
  replicates = object$replicates$parameters
  parameterNames = resolveUncertaintyParameters(names(replicates), parameters)
  showPoints = resolveUncertaintyShowPoints(showPoints, default = TRUE)

  plotUncertaintyReplicates(
    replicates = replicates,
    level = level,
    parameters = parameterNames,
    interpretation = "Rubin Bayesian Bootstrap weighted-fit uncertainty",
    showPoints = showPoints,
    xlab = xlab,
    ylab = ylab,
    main = main,
    ...
  )
}

#' Validate uncertainty probability levels.
#'
#' @param level Probability content levels.
#' @return Sorted unique levels.
#' @keywords internal
#' @noRd
validateUncertaintyLevels = function(level) {
  if (!is.numeric(level) || length(level) < 1L || any(!is.finite(level)) ||
      any(level <= 0 | level >= 1)) {
    stop("level must contain finite probability values strictly between 0 and 1", call. = FALSE)
  }
  sort(unique(as.numeric(level)))
}

#' Validate uncertainty plotting grid size.
#'
#' @param nGrid Grid size.
#' @return Validated integer grid size.
#' @keywords internal
#' @noRd
validateUncertaintyGridSize = function(nGrid) {
  if (!is.numeric(nGrid) || length(nGrid) != 1L || !is.finite(nGrid) ||
      nGrid < 51L || nGrid != floor(nGrid)) {
    stop("nGrid must be one integer greater than or equal to 51", call. = FALSE)
  }
  as.integer(nGrid)
}

#' Resolve whether sample realizations should be shown.
#'
#' @param showPoints User-supplied logical value or `NULL`.
#' @param default Method-specific default.
#' @return One logical value.
#' @keywords internal
#' @noRd
resolveUncertaintyShowPoints = function(showPoints, default) {
  if (is.null(showPoints)) {
    return(isTRUE(default))
  }
  if (!is.logical(showPoints) || length(showPoints) != 1L || is.na(showPoints)) {
    stop("showPoints must be TRUE, FALSE, or NULL", call. = FALSE)
  }
  showPoints
}

#' Resolve one or two parameter names for an uncertainty display.
#'
#' @param available Available parameter names.
#' @param requested Optional requested parameter names.
#' @return Character vector of one or two parameter names.
#' @keywords internal
#' @noRd
resolveUncertaintyParameters = function(available, requested = NULL) {
  available = as.character(available)
  if (length(available) < 1L) {
    stop("no parameter uncertainty is available to plot", call. = FALSE)
  }

  if (is.null(requested)) {
    if (length(available) > 2L) {
      stop(
        "models with more than two parameters require parameters = c(...) selecting one or two parameters",
        call. = FALSE
      )
    }
    return(available)
  }

  if (!is.character(requested) || length(requested) < 1L || length(requested) > 2L ||
      anyDuplicated(requested) || !all(requested %in% available)) {
    stop(
      "parameters must name one or two available model parameters: ",
      paste(available, collapse = ", "),
      call. = FALSE
    )
  }
  requested
}

#' Plot uncertainty from stored parameter replicates.
#'
#' @keywords internal
#' @noRd
plotUncertaintyReplicates = function(replicates,
                                      level,
                                      parameters,
                                      interpretation,
                                      showPoints,
                                      xlab,
                                      ylab,
                                      main,
                                      positiveKde = FALSE,
                                      ...) {
  level = validateUncertaintyLevels(level)
  values = as.data.frame(replicates)[, parameters, drop = FALSE]
  values = values[complete.cases(values), , drop = FALSE]
  if (nrow(values) < 2L) {
    stop("at least two successful parameter replicates are required", call. = FALSE)
  }

  if (length(parameters) == 1L) {
    parameter = parameters[[1L]]
    densityEstimate = density(values[[parameter]])
    intervals = t(vapply(level, function(probabilityLevel) {
      alpha = (1 - probabilityLevel) / 2
      quantile(
        values[[parameter]],
        probs = c(alpha, 1 - alpha),
        names = FALSE,
        type = 7
      )
    }, numeric(2L)))
    colnames(intervals) = c("lower", "upper")
    rownames(intervals) = as.character(level)

    if (is.null(xlab)) {
      xlab = parameter
    }
    if (is.null(ylab)) {
      ylab = "Density"
    }
    if (is.null(main)) {
      main = interpretation
    }

    plot(
      densityEstimate$x,
      densityEstimate$y,
      type = "l",
      xlab = xlab,
      ylab = ylab,
      main = main,
      ...
    )
    for (intervalIndex in seq_len(nrow(intervals))) {
      abline(v = intervals[intervalIndex, ], lty = intervalIndex + 1L)
    }

    return(invisible(list(
      method = interpretation,
      dimension = 1L,
      parameters = parameters,
      level = level,
      density = data.frame(x = densityEstimate$x, density = densityEstimate$y),
      intervals = intervals,
      replicates = values
    )))
  }

  region = makeSampleUncertaintyRegion(
    values = values,
    parameters = parameters,
    level = level,
    positive = positiveKde
  )
  if (is.null(xlab)) {
    xlab = parameters[[1L]]
  }
  if (is.null(ylab)) {
    ylab = parameters[[2L]]
  }
  if (is.null(main)) {
    main = interpretation
  }

  plot(
    range(values[[parameters[[1L]]]]),
    range(values[[parameters[[2L]]]]),
    type = "n",
    xlab = xlab,
    ylab = ylab,
    main = main,
    ...
  )
  if (showPoints) {
    points(values[[parameters[[1L]]]], values[[parameters[[2L]]]], pch = 1)
  }
  drawUncertaintyContours(region$contours, parameters)

  invisible(list(
    method = interpretation,
    dimension = 2L,
    parameters = parameters,
    level = level,
    density = region$density,
    contours = region$contours,
    containment = region$containment,
    replicates = values
  ))
}

#' Construct KDE probability-content contours from parameter replicates.
#'
#' @keywords internal
#' @noRd
makeSampleUncertaintyRegion = function(values,
                                        parameters,
                                        level,
                                        weights = NULL,
                                        positive = FALSE) {
  matrixValues = as.matrix(values[, parameters, drop = FALSE])
  bandwidth = Hscv(matrixValues)
  densityArguments = list(x = matrixValues, H = bandwidth)
  if (!is.null(weights)) {
    weights = as.numeric(weights)
    if (length(weights) != nrow(matrixValues) || any(!is.finite(weights)) ||
        any(weights < 0) || sum(weights) <= 0) {
      stop("sample weights must be finite, non-negative, and match parameter replicates", call. = FALSE)
    }
    densityArguments$w = weights / sum(weights) * nrow(matrixValues)
  }
  if (!is.logical(positive) || length(positive) != 1L || is.na(positive)) {
    stop("positive must be TRUE or FALSE", call. = FALSE)
  }
  if (positive) {
    densityArguments$positive = TRUE
  }
  densityEstimate = do.call(kde, densityArguments)
  contourHeights = contourLevels(
    densityEstimate,
    cont = sort(100 * level),
    approx = TRUE
  )
  contours = contourLines(
    x = densityEstimate$eval.points[[1L]],
    y = densityEstimate$eval.points[[2L]],
    z = densityEstimate$estimate,
    levels = contourHeights
  )
  contours = labelUncertaintyContours(contours, contourHeights, level, parameters)
  containment = uncertaintyContourContainment(
    values = matrixValues,
    contours = contours,
    weights = weights
  )

  list(
    density = densityEstimate,
    contours = contours,
    containment = containment
  )
}

#' Calculate empirical containment of sample-based uncertainty contours.
#'
#' @param values Matrix of two-dimensional parameter realizations.
#' @param contours Labeled contour coordinates.
#' @param weights Optional non-negative realization weights.
#' @return A data frame with nominal level and empirical containment.
#' @keywords internal
#' @noRd
uncertaintyContourContainment = function(values, contours, weights = NULL) {
  values = as.matrix(values)
  if (ncol(values) != 2L || nrow(values) < 1L) {
    stop("values must contain at least one two-dimensional realization", call. = FALSE)
  }
  if (length(contours) == 0L) {
    return(data.frame(level = numeric(0), containment = numeric(0)))
  }

  if (is.null(weights)) {
    weights = rep(1, nrow(values))
  } else {
    weights = as.numeric(weights)
    if (length(weights) != nrow(values) || any(!is.finite(weights)) ||
        any(weights < 0) || sum(weights) <= 0) {
      stop("sample weights must be finite, non-negative, and match parameter realizations", call. = FALSE)
    }
  }
  weights = weights / sum(weights)

  levels = sort(unique(vapply(contours, `[[`, numeric(1L), "level")))
  containment = vapply(levels, function(probabilityLevel) {
    levelContours = contours[vapply(contours, function(contour) {
      isTRUE(all.equal(contour$level, probabilityLevel))
    }, logical(1L))]
    inside = rep(FALSE, nrow(values))
    for (contour in levelContours) {
      inside = inside | pointsInsidePolygon(
        x = values[, 1L],
        y = values[, 2L],
        polygonX = contour$x,
        polygonY = contour$y
      )
    }
    sum(weights[inside])
  }, numeric(1L))

  data.frame(level = levels, containment = containment)
}

#' Determine whether points lie inside a polygon.
#'
#' Uses the standard ray-crossing rule and treats points on polygon edges as
#' inside. This keeps contour-containment diagnostics self-contained without
#' adding a spatial package dependency.
#'
#' @param x Point x coordinates.
#' @param y Point y coordinates.
#' @param polygonX Polygon x coordinates.
#' @param polygonY Polygon y coordinates.
#' @return Logical vector indicating polygon membership.
#' @keywords internal
#' @noRd
pointsInsidePolygon = function(x, y, polygonX, polygonY) {
  x = as.numeric(x)
  y = as.numeric(y)
  polygonX = as.numeric(polygonX)
  polygonY = as.numeric(polygonY)
  if (length(x) != length(y)) {
    stop("x and y must have the same length", call. = FALSE)
  }
  if (length(polygonX) != length(polygonY) || length(polygonX) < 3L) {
    stop("polygon coordinates must have equal lengths of at least three", call. = FALSE)
  }

  if (!isTRUE(all.equal(polygonX[[1L]], polygonX[[length(polygonX)]])) ||
      !isTRUE(all.equal(polygonY[[1L]], polygonY[[length(polygonY)]]))) {
    polygonX = c(polygonX, polygonX[[1L]])
    polygonY = c(polygonY, polygonY[[1L]])
  }

  inside = rep(FALSE, length(x))
  onBoundary = rep(FALSE, length(x))
  tolerance = sqrt(.Machine$double.eps)

  for (edgeIndex in seq_len(length(polygonX) - 1L)) {
    x1 = polygonX[[edgeIndex]]
    y1 = polygonY[[edgeIndex]]
    x2 = polygonX[[edgeIndex + 1L]]
    y2 = polygonY[[edgeIndex + 1L]]

    crossProduct = (x - x1) * (y2 - y1) - (y - y1) * (x2 - x1)
    edgeScale = max(abs(c(x1, y1, x2, y2)), 1)
    withinX = x >= min(x1, x2) - tolerance * edgeScale &
      x <= max(x1, x2) + tolerance * edgeScale
    withinY = y >= min(y1, y2) - tolerance * edgeScale &
      y <= max(y1, y2) + tolerance * edgeScale
    onBoundary = onBoundary |
      (abs(crossProduct) <= tolerance * edgeScale & withinX & withinY)

    crosses = (y1 > y) != (y2 > y)
    crossingX = rep(Inf, length(x))
    if (y2 != y1 && any(crosses)) {
      crossingX[crosses] = x1 +
        (y[crosses] - y1) * (x2 - x1) / (y2 - y1)
    }
    inside = xor(inside, crosses & x < crossingX)
  }

  inside | onBoundary
}

#' Label contour-line output with probability content and parameter names.
#'
#' @keywords internal
#' @noRd
labelUncertaintyContours = function(contours, heights, level, parameters) {
  if (length(contours) == 0L) {
    return(list())
  }

  lapply(contours, function(contour) {
    levelIndex = which.min(abs(heights - contour$level))
    list(
      level = level[[levelIndex]],
      densityLevel = contour$level,
      x = contour$x,
      y = contour$y,
      parameters = parameters
    )
  })
}

#' Draw stored uncertainty contour coordinates.
#'
#' @keywords internal
#' @noRd
drawUncertaintyContours = function(contours, parameters) {
  if (length(contours) == 0L) {
    return(invisible(NULL))
  }
  uniqueLevels = sort(unique(vapply(contours, `[[`, numeric(1L), "level")))
  for (contour in contours) {
    lineType = match(contour$level, uniqueLevels)
    lines(contour$x, contour$y, lty = lineType, lwd = 2)
  }
  invisible(NULL)
}

#' Plot maximum-likelihood parameter uncertainty.
#'
#' @keywords internal
#' @noRd
plotUncertaintyMle = function(object,
                               model,
                               level,
                               parameters,
                               nGrid,
                               xlab,
                               ylab,
                               main,
                               ...) {
  level = validateUncertaintyLevels(level)

  if (length(parameters) == 2L) {
    if (!identical(object$model, "ziz") ||
        !identical(parameters, c("pi", "shape"))) {
      stop(
        "two-parameter MLE uncertainty plotting currently uses the established ZIZ profile-likelihood region",
        call. = FALSE
      )
    }
    regions = profileLikelihoodZIZ(object, level = level, silent = TRUE)
    if (is.null(xlab)) {
      xlab = "pi"
    }
    if (is.null(ylab)) {
      ylab = "shape"
    }
    if (is.null(main)) {
      main = "Profile-likelihood confidence region"
    }

    allPi = unlist(lapply(regions, `[[`, "pi"))
    allShape = unlist(lapply(regions, `[[`, "shape"))
    plot(
      range(allPi),
      range(allShape),
      type = "n",
      xlab = xlab,
      ylab = ylab,
      main = main,
      ...
    )
    for (regionIndex in seq_along(regions)) {
      lines(regions[[regionIndex]]$pi, regions[[regionIndex]]$shape, lty = regionIndex, lwd = 2)
    }

    return(invisible(list(
      method = "Profile-likelihood confidence",
      dimension = 2L,
      parameters = parameters,
      level = level,
      contours = regions
    )))
  }

  profile = makeOneDimensionalLikelihoodProfile(
    object = object,
    model = model,
    parameter = parameters[[1L]],
    level = level,
    nGrid = nGrid
  )
  if (is.null(xlab)) {
    xlab = parameters[[1L]]
  }
  if (is.null(ylab)) {
    ylab = "Relative likelihood"
  }
  if (is.null(main)) {
    main = "Profile-likelihood confidence interval"
  }

  plot(
    profile$x,
    profile$relativeLikelihood,
    type = "l",
    xlab = xlab,
    ylab = ylab,
    main = main,
    ...
  )
  for (intervalIndex in seq_len(nrow(profile$intervals))) {
    abline(v = profile$intervals[intervalIndex, ], lty = intervalIndex + 1L)
  }

  invisible(list(
    method = "Profile-likelihood confidence",
    dimension = 1L,
    parameters = parameters,
    level = level,
    density = data.frame(
      x = profile$x,
      relativeLikelihood = profile$relativeLikelihood
    ),
    intervals = profile$intervals
  ))
}

#' Build a one-dimensional relative-likelihood profile.
#'
#' @keywords internal
#' @noRd
makeOneDimensionalLikelihoodProfile = function(object,
                                                model,
                                                parameter,
                                                level,
                                                nGrid) {
  estimate = fitModelParameters(object, model)[[parameter]]
  spread = getMleParameterSpread(object, parameter)
  control = modelMleControl(model, object$psData)
  lower = unname(control$lower[[parameter]])
  upper = unname(control$upper[[parameter]])

  if (!is.finite(spread) || spread <= 0) {
    spread = max(abs(estimate) / 4, 0.1)
  }
  gridLower = if (is.finite(lower)) {
    max(lower + sqrt(.Machine$double.eps), estimate - 6 * spread)
  } else {
    estimate - 6 * spread
  }
  gridUpper = if (is.finite(upper)) {
    min(upper - sqrt(.Machine$double.eps), estimate + 6 * spread)
  } else {
    estimate + 6 * spread
  }
  if (!is.finite(gridLower) || !is.finite(gridUpper) || gridLower >= gridUpper) {
    stop("unable to construct a useful likelihood plotting range", call. = FALSE)
  }

  x = seq(gridLower, gridUpper, length.out = as.integer(nGrid))
  logLikelihood = vapply(x, function(parameterValue) {
    parameterValues = c(parameterValue)
    names(parameterValues) = parameter
    modelLogLikelihood(model, parameters = parameterValues, data = object$psData)
  }, numeric(1L))
  relativeLikelihood = exp(logLikelihood - max(logLikelihood, na.rm = TRUE))

  intervals = t(vapply(level, function(probabilityLevel) {
    threshold = exp(-0.5 * qchisq(probabilityLevel, df = 1L))
    included = x[relativeLikelihood >= threshold]
    if (length(included) == 0L) {
      return(c(NA_real_, NA_real_))
    }
    range(included)
  }, numeric(2L)))
  colnames(intervals) = c("lower", "upper")
  rownames(intervals) = as.character(level)

  list(
    x = x,
    relativeLikelihood = relativeLikelihood,
    intervals = intervals
  )
}

#' Extract an MLE parameter standard error for plotting-range selection.
#'
#' @keywords internal
#' @noRd
getMleParameterSpread = function(object, parameter) {
  if (identical(parameter, "shape") && !is.null(object$var.shape)) {
    return(sqrt(as.numeric(object$var.shape)))
  }
  if (identical(parameter, "pi") && !is.null(object$var.pi)) {
    return(sqrt(as.numeric(object$var.pi)))
  }
  if (is.matrix(object$var.cov) && parameter %in% rownames(object$var.cov)) {
    return(sqrt(as.numeric(object$var.cov[parameter, parameter])))
  }
  NA_real_
}

#' Plot parametric Bayesian parameter uncertainty.
#'
#' @keywords internal
#' @noRd
plotUncertaintyBayes = function(object,
                                 level,
                                 parameters,
                                 showPoints,
                                 nGrid,
                                 xlab,
                                 ylab,
                                 main,
                                 ...) {
  if (!inherits(object$posterior, "psPosterior")) {
    stop("Bayesian psFit object does not contain a psPosterior representation", call. = FALSE)
  }

  if (length(parameters) == 1L) {
    posteriorData = getPosteriorPlotData(
      object = object,
      parameter = parameters[[1L]],
      level = max(level),
      nGrid = nGrid
    )
    intervals = t(vapply(level, function(probabilityLevel) {
      posteriorLevelData = getPosteriorPlotData(
        object = object,
        parameter = parameters[[1L]],
        level = probabilityLevel,
        nGrid = nGrid
      )
      posteriorLevelData$interval
    }, numeric(2L)))
    colnames(intervals) = c("lower", "upper")
    rownames(intervals) = as.character(level)

    if (is.null(xlab)) {
      xlab = parameters[[1L]]
    }
    if (is.null(ylab)) {
      ylab = "Posterior density"
    }
    if (is.null(main)) {
      main = "Parametric Bayesian posterior uncertainty"
    }
    plot(
      posteriorData$x,
      posteriorData$density,
      type = "l",
      xlab = xlab,
      ylab = ylab,
      main = main,
      ...
    )
    for (intervalIndex in seq_len(nrow(intervals))) {
      abline(v = intervals[intervalIndex, ], lty = intervalIndex + 1L)
    }

    return(invisible(list(
      method = "Parametric Bayesian credible",
      dimension = 1L,
      parameters = parameters,
      level = level,
      density = data.frame(x = posteriorData$x, density = posteriorData$density),
      intervals = intervals
    )))
  }

  representation = object$posterior$representation
  chain = representation$value$chain
  if (!is.null(chain)) {
    return(plotUncertaintyReplicates(
      replicates = as.data.frame(chain),
      level = level,
      parameters = parameters,
      interpretation = "Parametric Bayesian posterior credible",
      showPoints = showPoints,
      xlab = xlab,
      ylab = ylab,
      main = main,
      positiveKde = TRUE,
      ...
    ))
  }

  approximation = representation$value$approximation
  if (!is.null(approximation$samples) &&
      is.data.frame(approximation$samples) &&
      all(c(parameters, "weight") %in% names(approximation$samples))) {
    values = approximation$samples[, parameters, drop = FALSE]
    region = makeSampleUncertaintyRegion(
      values = values,
      parameters = parameters,
      level = level,
      weights = approximation$samples$weight,
      positive = TRUE
    )
    if (is.null(xlab)) {
      xlab = parameters[[1L]]
    }
    if (is.null(ylab)) {
      ylab = parameters[[2L]]
    }
    if (is.null(main)) {
      main = "Importance-sampling posterior credible region"
    }

    plot(
      range(values[[parameters[[1L]]]]),
      range(values[[parameters[[2L]]]]),
      type = "n",
      xlab = xlab,
      ylab = ylab,
      main = main,
      ...
    )
    if (showPoints) {
      points(values[[parameters[[1L]]]], values[[parameters[[2L]]]], pch = 1)
    }
    drawUncertaintyContours(region$contours, parameters)

    return(invisible(list(
      method = "Parametric Bayesian credible",
      dimension = 2L,
      parameters = parameters,
      level = level,
      density = region$density,
      contours = region$contours,
      replicates = values,
      weights = approximation$samples$weight
    )))
  }

  if (!is.null(approximation$mode) && !is.null(approximation$varCov)) {
    contours = makeGaussianUncertaintyContours(
      centre = approximation$mode[parameters],
      covariance = approximation$varCov[parameters, parameters, drop = FALSE],
      level = level,
      parameters = parameters
    )
    if (is.null(xlab)) {
      xlab = parameters[[1L]]
    }
    if (is.null(ylab)) {
      ylab = parameters[[2L]]
    }
    if (is.null(main)) {
      main = "Laplace posterior credible region"
    }

    allX = unlist(lapply(contours, `[[`, "x"))
    allY = unlist(lapply(contours, `[[`, "y"))
    plot(
      range(allX),
      range(allY),
      type = "n",
      xlab = xlab,
      ylab = ylab,
      main = main,
      ...
    )
    drawUncertaintyContours(contours, parameters)

    return(invisible(list(
      method = "Parametric Bayesian credible",
      dimension = 2L,
      parameters = parameters,
      level = level,
      contours = contours,
      approximation = "laplace"
    )))
  }

  grid = representation$value$grid
  if (!is.null(grid) && !is.null(grid$parameters) && !is.null(grid$weights)) {
    region = makeGridUncertaintyRegion(grid, parameters, level)
    if (is.null(xlab)) {
      xlab = parameters[[1L]]
    }
    if (is.null(ylab)) {
      ylab = parameters[[2L]]
    }
    if (is.null(main)) {
      main = "Numerical posterior credible region"
    }

    plot(
      range(grid$parameters[[parameters[[1L]]]]),
      range(grid$parameters[[parameters[[2L]]]]),
      type = "n",
      xlab = xlab,
      ylab = ylab,
      main = main,
      ...
    )
    drawUncertaintyContours(region$contours, parameters)

    return(invisible(list(
      method = "Parametric Bayesian credible",
      dimension = 2L,
      parameters = parameters,
      level = level,
      cumulativeMass = region$cumulativeMass,
      contours = region$contours,
      grid = grid
    )))
  }

  stop(
    "two-dimensional plotUncertainty() could not extract a supported posterior uncertainty representation",
    call. = FALSE
  )
}

#' Construct Gaussian probability contours for a Laplace approximation.
#'
#' @keywords internal
#' @noRd
makeGaussianUncertaintyContours = function(centre, covariance, level, parameters) {
  centre = as.numeric(centre)
  covariance = as.matrix(covariance)
  if (length(centre) != 2L || !identical(dim(covariance), c(2L, 2L)) ||
      any(!is.finite(centre)) || any(!is.finite(covariance))) {
    stop("Laplace uncertainty requires a finite two-parameter centre and covariance", call. = FALSE)
  }

  decomposition = eigen(covariance, symmetric = TRUE)
  if (any(!is.finite(decomposition$values)) || any(decomposition$values <= 0)) {
    stop("Laplace covariance must be positive definite", call. = FALSE)
  }
  transform = decomposition$vectors %*% diag(sqrt(decomposition$values), nrow = 2L)
  angles = seq(0, 2 * pi, length.out = 241L)
  unitCircle = rbind(cos(angles), sin(angles))

  lapply(level, function(probabilityLevel) {
    radius = sqrt(qchisq(probabilityLevel, df = 2L))
    coordinates = centre + transform %*% (radius * unitCircle)
    list(
      level = probabilityLevel,
      densityLevel = NA_real_,
      x = coordinates[1L, ],
      y = coordinates[2L, ],
      parameters = parameters
    )
  })
}

#' Construct probability-content contours from a stored numerical posterior grid.
#'
#' The stored grid masses are sorted from largest to smallest and cumulatively
#' summed. The density/mass threshold at which each requested probability level
#' is reached defines the corresponding highest-density contour.
#'
#' @keywords internal
#' @noRd
makeGridUncertaintyRegion = function(grid, parameters, level) {
  parameterGrid = as.data.frame(grid$parameters)
  mass = if (!is.null(grid$mass)) {
    as.numeric(grid$mass)
  } else {
    as.numeric(grid$weights)
  }
  densityValues = if (!is.null(grid$density)) {
    as.numeric(grid$density)
  } else {
    mass
  }

  if (nrow(parameterGrid) != length(mass) ||
      nrow(parameterGrid) != length(densityValues) ||
      any(!is.finite(mass)) || any(mass < 0) || sum(mass) <= 0 ||
      any(!is.finite(densityValues)) || any(densityValues < 0)) {
    stop("numerical posterior grid is incomplete", call. = FALSE)
  }
  mass = mass / sum(mass)

  ordering = order(densityValues, decreasing = TRUE)
  cumulativeMass = cumsum(mass[ordering])
  thresholds = vapply(level, function(probabilityLevel) {
    index = which(cumulativeMass >= probabilityLevel)[1L]
    densityValues[ordering[[index]]]
  }, numeric(1L))

  xValues = sort(unique(parameterGrid[[parameters[[1L]]]]))
  yValues = sort(unique(parameterGrid[[parameters[[2L]]]]))
  z = matrix(
    NA_real_,
    nrow = length(xValues),
    ncol = length(yValues)
  )
  xIndex = match(parameterGrid[[parameters[[1L]]]], xValues)
  yIndex = match(parameterGrid[[parameters[[2L]]]], yValues)
  z[cbind(xIndex, yIndex)] = densityValues

  rawContours = contourLines(
    x = xValues,
    y = yValues,
    z = z,
    levels = sort(unique(thresholds))
  )
  contours = lapply(rawContours, function(contour) {
    thresholdIndex = which.min(abs(thresholds - contour$level))
    list(
      level = level[[thresholdIndex]],
      densityLevel = contour$level,
      x = contour$x,
      y = contour$y,
      parameters = parameters
    )
  })

  list(
    contours = contours,
    thresholds = thresholds,
    cumulativeMass = data.frame(
      rank = seq_along(ordering),
      density = densityValues[ordering],
      mass = mass[ordering],
      cumulative = cumulativeMass
    )
  )
}
