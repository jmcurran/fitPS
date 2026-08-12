#' S3 predict method for an object of class \code{psFit}
#'
#' @param object an object of class \code{psFit}, usually from \code{\link{fitDist}}
#' or \code{\link{fitZIDist}}.
#' @param newdata an optional vector of integers at which to calculate
#'   \eqn{\Pr(X = x)}{Pr(X = x)}.
#' @param type The probability definition. \code{"plugIn"} evaluates the
#'   probability model at fitted parameter values. \code{"posteriorMean"}
#'   returns posterior mean probabilities for Bayesian fits.
#'   \code{"bootstrapMean"} returns bootstrap mean probabilities when a
#'   bootstrap distribution has been attached to the fit.
#' @param interval Interval type. Existing \code{"prof"} and \code{"wald"}
#'   intervals are retained for plug-in zeta predictions. \code{"credible"}
#'   returns stored equal-tailed credible intervals for posterior mean
#'   predictions. \code{"percentile"} returns stored percentile bootstrap
#'   confidence intervals for bootstrap mean predictions.
#' @param level the interval level. For posterior and bootstrap predictions,
#'   the requested level must match the level stored in the corresponding
#'   uncertainty object. If omitted, the stored level is used.
#' @param ... other arguments passed to \code{predict}---not used.
#'
#' @return either a named vector of probabilities, or a \code{data.frame}
#'   with columns \code{predicted}, \code{lower}, and \code{upper} and row
#'   names showing the requested P or S terms.
#' @examples
#' data(Psurveys)
#' roux = Psurveys$roux
#' fit = fitDist(roux)
#' predict(fit, interval = "prof")
#' @importFrom stats predict qnorm
#' @export
predict.psFit = function(object,
                         newdata,
                         type = c("plugIn", "posteriorMean", "bootstrapMean"),
                         interval = c(
                           "none", "prof", "wald", "credible", "percentile"
                         ),
                         level = 0.95,
                         ...) {
  type = match.arg(type)
  interval = match.arg(interval)
  levelMissing = missing(level)

  if (missing(newdata)) {
    newdata = as.numeric(
      gsub("^(P|S)([0-9]+)$", "\\2", names(object$fitted))
    )
  }

  validatePredictionIndices(object, newdata)

  if (identical(type, "posteriorMean")) {
    if (!inherits(object$posterior, "psPosterior")) {
      stop(
        "posteriorMean predictions require a Bayesian psFit object with ",
        "posterior summaries"
      )
    }
    if (!interval %in% c("none", "credible")) {
      stop(
        "posteriorMean predictions support interval = 'none' or 'credible'"
      )
    }

    if (levelMissing) {
      level = object$posterior$level
    }
    validateStoredPredictionLevel(level, object$posterior$level, "posterior")

    summaries = selectPredictionSummaries(
      object$posterior$probabilities,
      newdata,
      object$psData$type
    )
    return(formatInferencePredictions(
      summaries,
      interval = interval,
      intervalName = "credible"
    ))
  }

  if (identical(type, "bootstrapMean")) {
    if (!inherits(object$bootstrap, "psBootstrap")) {
      stop(
        "bootstrapMean predictions require a psFit object with an attached ",
        "bootstrap; run bootstrapFit() first"
      )
    }
    if (!interval %in% c("none", "percentile")) {
      stop(
        "bootstrapMean predictions support interval = 'none' or 'percentile'"
      )
    }

    if (levelMissing) {
      level = object$bootstrap$level
    }
    validateStoredPredictionLevel(level, object$bootstrap$level, "bootstrap")

    summaries = selectPredictionSummaries(
      object$bootstrap$probabilities,
      newdata,
      object$psData$type
    )
    return(formatInferencePredictions(
      summaries,
      interval = interval,
      intervalName = "percentile"
    ))
  }

  if (interval %in% c("credible", "percentile")) {
    stop(
      "plugIn predictions do not support credible or percentile intervals; ",
      "use type = 'posteriorMean' or type = 'bootstrapMean'"
    )
  }

  predicted = calculatePlugInPredictions(object, newdata)

  if (interval %in% c("prof", "wald")) {
    if (object$model == "zeta") {
      if (level <= 0.75 || level >= 1) {
        stop("Level should be in the interval [0.75, 1)")
      }

      zstar = qnorm((1 - level) * 0.5, lower.tail = FALSE)
      offset = ifelse(object$psData$type == "P", 1, 0)
      lower = dzetaStandard(
        newdata + offset,
        shape = object$shape - zstar * sqrt(object$var.shape)
      )
      upper = dzetaStandard(
        newdata + offset,
        shape = object$shape + zstar * sqrt(object$var.shape)
      )

      result = data.frame(
        predicted = predicted,
        lower = lower,
        upper = upper
      )
      rownames(result) = paste0(object$psData$type, newdata)
      return(result)
    }

    names(predicted) = paste0(object$psData$type, newdata)
    return(predicted)
  }

  names(predicted) = paste0(object$psData$type, newdata)
  predicted
}

#' Validate prediction indices
#'
#' @param object A `psFit` object.
#' @param newdata Numeric prediction indices.
#' @return `newdata`, invisibly.
#' @noRd
validatePredictionIndices = function(object, newdata) {
  if (!is.numeric(newdata) || length(newdata) < 1L || any(!is.finite(newdata))) {
    stop("newdata should contain finite numeric values")
  }

  if (any(abs(newdata - floor(newdata)) > .Machine$double.eps^0.5)) {
    stop("newdata should only contain integers")
  }

  if (object$psData$type == "P" && any(newdata < 0)) {
    stop("Can only make predictions for P probabilities for values of n >= 0")
  }

  if (object$psData$type == "S" && any(newdata <= 0)) {
    stop("Can only make predictions for size probabilities for values of n >= 1")
  }

  invisible(newdata)
}

#' Calculate plug-in predictions
#'
#' @param object A `psFit` object.
#' @param newdata Integer prediction indices.
#' @return A numeric vector of plug-in probabilities.
#' @noRd
calculatePlugInPredictions = function(object, newdata) {
  if (object$model %in% c("zeta", "ziz")) {
    probabilityFunction = probfun(object)
    return(unname(probabilityFunction(newdata)))
  }

  stop("This method is not currently implemented for the logarithmic distribution")
}

#' Select stored uncertainty summaries for prediction
#'
#' @param probabilities A posterior or bootstrap probability summary data frame.
#' @param newdata Integer prediction indices.
#' @param surveyType Either `"P"` or `"S"`.
#' @return The selected probability summary rows.
#' @noRd
selectPredictionSummaries = function(probabilities, newdata, surveyType) {
  requestedTerms = paste0(surveyType, as.integer(newdata))
  selectedRows = match(requestedTerms, probabilities$term)

  if (anyNA(selectedRows)) {
    missingTerms = requestedTerms[is.na(selectedRows)]
    stop(
      "stored uncertainty summaries are not available for: ",
      paste(missingTerms, collapse = ", "),
      "; refit with enough terms before requesting these predictions"
    )
  }

  probabilities[selectedRows, , drop = FALSE]
}

#' Validate an interval level against stored summaries
#'
#' @param level Requested interval level.
#' @param storedLevel Interval level stored in the uncertainty object.
#' @param label Character label used in error messages.
#' @return `level`, invisibly.
#' @noRd
validateStoredPredictionLevel = function(level, storedLevel, label) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("level must be one finite number strictly between 0 and 1")
  }

  if (!isTRUE(all.equal(level, storedLevel))) {
    stop(
      label,
      " summaries are stored at level ",
      format(storedLevel),
      "; refit or recompute the uncertainty object to use level ",
      format(level)
    )
  }

  invisible(level)
}

#' Format posterior or bootstrap predictions
#'
#' @param summaries A probability summary data frame.
#' @param interval Requested interval type.
#' @param intervalName Interval type supported by these summaries.
#' @return A named numeric vector or prediction interval data frame.
#' @noRd
formatInferencePredictions = function(summaries,
                                       interval,
                                       intervalName) {
  estimates = summaries$estimate
  names(estimates) = summaries$term

  if (identical(interval, "none")) {
    return(estimates)
  }

  if (!identical(interval, intervalName)) {
    stop("invalid interval type for the requested prediction definition")
  }

  result = data.frame(
    predicted = summaries$estimate,
    lower = summaries$lower,
    upper = summaries$upper
  )
  rownames(result) = summaries$term
  result
}
