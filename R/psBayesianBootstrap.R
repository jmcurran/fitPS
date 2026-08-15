#' Construct an internal Rubin Bayesian Bootstrap object
#'
#' @param method Character description of the Bayesian Bootstrap method.
#' @param parameters Data frame of parameter summaries.
#' @param probabilities Data frame of fitted-probability summaries.
#' @param replicates List containing parameter and probability replicate data.
#' @param level Probability level used for equal-tail intervals.
#' @param diagnostics Replicate diagnostics, including failed-fit information.
#' @param seed Optional random-number seed used to generate the draws.
#' @param model Model identifier associated with the weighted refits.
#' @return An internal object of class `psBayesianBootstrap`.
#' @importFrom stats quantile rgamma sd
#' @keywords internal
#' @noRd
newPsBayesianBootstrap = function(method,
                                   parameters,
                                   probabilities,
                                   replicates,
                                   level,
                                   diagnostics,
                                   seed = NULL,
                                   model = NULL) {
  if (!is.character(method) || length(method) != 1L || !nzchar(method)) {
    stop("method must be one non-empty character value")
  }
  if (!is.data.frame(parameters) ||
      !all(c("parameter", "estimate", "sd", "lower", "upper", "level") %in%
        names(parameters))) {
    stop("parameters must contain Bayesian Bootstrap parameter summaries")
  }
  if (!is.data.frame(probabilities) ||
      !all(c("term", "estimate", "sd", "lower", "upper", "level") %in%
        names(probabilities))) {
    stop("probabilities must contain Bayesian Bootstrap probability summaries")
  }
  if (!is.list(replicates) ||
      !all(c("parameters", "probabilities") %in% names(replicates)) ||
      !is.data.frame(replicates$parameters) ||
      !is.data.frame(replicates$probabilities)) {
    stop("replicates must contain parameter and probability data frames")
  }

  result = list(
    method = method,
    parameters = parameters,
    probabilities = probabilities,
    replicates = replicates,
    level = level,
    diagnostics = diagnostics,
    seed = seed,
    model = model
  )
  class(result) = "psBayesianBootstrap"
  result
}

#' Validate a Bayesian Bootstrap replicate count
#'
#' @param B Number of requested Bayesian Bootstrap replicates.
#' @return The validated integer replicate count.
#' @keywords internal
#' @noRd
validateBayesianBootstrapB = function(B) {
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) ||
      B < 1 || B != floor(B)) {
    stop("B must be one positive integer")
  }
  as.integer(B)
}

#' Validate an optional Bayesian Bootstrap random-number seed
#'
#' @param seed Optional finite numeric seed.
#' @return The validated integer seed or `NULL`.
#' @keywords internal
#' @noRd
validateBayesianBootstrapSeed = function(seed) {
  if (is.null(seed)) {
    return(NULL)
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
    stop("seed must be NULL or one finite numeric value")
  }
  as.integer(seed)
}

#' Draw grouped Rubin Bayesian Bootstrap weights
#'
#' For an observed survey with occupied-category counts `n_j`, the sums of
#' individual `Dirichlet(1, ..., 1)` weights within categories follow
#' `Dirichlet(n_1, ..., n_k)`. This helper draws those grouped weights directly
#' and scales each row to the original observed sample size.
#'
#' @param x An observed `psData` object.
#' @param B Number of weight draws.
#' @param seed Optional random-number seed.
#' @param sampleSizeScale Logical; if `TRUE`, scale each draw to sum to the
#'   observed sample size. Otherwise each draw sums to one.
#' @return A numeric matrix with one Bayesian Bootstrap draw per row.
#' @keywords internal
#' @noRd
drawGroupedBayesianBootstrapWeights = function(x,
                                                 B,
                                                 seed = NULL,
                                                 sampleSizeScale = TRUE) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData")
  }
  B = validateBayesianBootstrapB(B)
  seed = validateBayesianBootstrapSeed(seed)
  if (!is.logical(sampleSizeScale) || length(sampleSizeScale) != 1L ||
      is.na(sampleSizeScale)) {
    stop("sampleSizeScale must be TRUE or FALSE")
  }

  counts = x$data$rn
  if (!is.numeric(counts) || any(!is.finite(counts)) || any(counts <= 0)) {
    stop("psData occupied-category counts must be finite positive values")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  categoryCount = length(counts)
  gammaDraws = rgamma(
    B * categoryCount,
    shape = rep(counts, times = B),
    rate = 1
  )
  weightMatrix = matrix(
    gammaDraws,
    nrow = B,
    ncol = categoryCount,
    byrow = TRUE
  )
  weightMatrix = weightMatrix / rowSums(weightMatrix)

  if (sampleSizeScale) {
    weightMatrix = weightMatrix * sum(counts)
  }

  colnames(weightMatrix) = as.character(x$data$n)
  weightMatrix
}

#' Summarise Bayesian Bootstrap replicate columns
#'
#' @param replicates Data frame containing numeric replicate columns.
#' @param nameColumn Name of the output column that identifies each quantity.
#' @param level Probability level used for equal-tail intervals.
#' @return A data frame of means, standard deviations, and equal-tail bounds.
#' @keywords internal
#' @noRd
summariseBayesianBootstrapReplicates = function(replicates,
                                                 nameColumn,
                                                 level) {
  validateBootstrapLevel(level)

  summaries = lapply(names(replicates), function(columnName) {
    values = replicates[[columnName]]
    values = values[is.finite(values)]

    if (length(values) == 0L) {
      estimate = NA_real_
      standardDeviation = NA_real_
      bounds = c(NA_real_, NA_real_)
    } else {
      estimate = mean(values)
      standardDeviation = if (length(values) > 1L) {
        sd(values)
      } else {
        NA_real_
      }
      alpha = (1 - level) / 2
      bounds = unname(quantile(
        values,
        probs = c(alpha, 1 - alpha),
        names = FALSE,
        type = 7
      ))
    }

    result = data.frame(
      name = columnName,
      estimate = estimate,
      sd = standardDeviation,
      lower = bounds[1L],
      upper = bounds[2L],
      level = level,
      stringsAsFactors = FALSE
    )
    names(result)[1L] = nameColumn
    result
  })

  do.call(rbind, summaries)
}

#' Fit one Rubin Bayesian Bootstrap weighted replicate
#'
#' @param x Observed `psData` object.
#' @param model Built-in fitPS model descriptor.
#' @param weights Positive grouped Bayesian Bootstrap weights.
#' @param nterms Number of model-implied P/S probabilities to retain.
#' @return A list containing fitted parameters and probabilities, or a failure
#'   message if the weighted MLE cannot be obtained.
#' @keywords internal
#' @noRd
fitBayesianBootstrapReplicate = function(x, model, weights, nterms) {
  tryCatch(
    {
      weightedFit = fitWeightedModel(
        x = x,
        model = model,
        weights = weights,
        nterms = nterms
      )
      list(
        parameters = weightedFit$parameters,
        probabilities = weightedFit$fitted,
        error = NULL
      )
    },
    error = function(error) {
      list(
        parameters = NULL,
        probabilities = NULL,
        error = conditionMessage(error)
      )
    }
  )
}

#' Generate model-based Rubin Bayesian Bootstrap replicates
#'
#' The empirical support remains fixed. Each replicate draws grouped
#' `Dirichlet(n_1, ..., n_k)` weights, refits the selected model by weighted
#' MLE, and stores model-implied P/S probabilities. Failed fits remain present
#' as missing replicate rows and are described in diagnostics.
#'
#' @param x Observed `psData` object.
#' @param model Built-in fitPS model descriptor.
#' @param B Number of Bayesian Bootstrap replicates.
#' @param seed Optional random-number seed.
#' @param nterms Number of model-implied P/S probabilities to retain.
#' @param level Probability level used for equal-tail summaries.
#' @return An internal `psBayesianBootstrap` object.
#' @keywords internal
#' @noRd
bayesianBootstrapModel = function(x,
                                   model,
                                   B = 2000,
                                   seed = NULL,
                                   nterms = 10,
                                   level = 0.95) {
  if (!is(x, "psData")) {
    stop("x must be an object of class psData")
  }
  if (!inherits(model, "psModel")) {
    stop("model must inherit from psModel")
  }
  B = validateBayesianBootstrapB(B)
  seed = validateBayesianBootstrapSeed(seed)
  validateBootstrapLevel(level)
  if (!is.numeric(nterms) || length(nterms) != 1L || !is.finite(nterms) ||
      nterms < 1 || nterms != floor(nterms)) {
    stop("nterms must be one positive integer")
  }
  nterms = as.integer(nterms)

  parameterNames = modelParameterNames(model)
  probabilityIndices = posteriorProbabilityIndices(x$type, nterms)
  probabilityTemplate = modelProbabilities(
    model = model,
    parameters = as.list(weightedMleControl(model)$start),
    n = probabilityIndices,
    type = x$type
  )
  probabilityNames = colnames(probabilityTemplate)

  weightDraws = drawGroupedBayesianBootstrapWeights(
    x = x,
    B = B,
    seed = seed,
    sampleSizeScale = TRUE
  )

  parameterMatrix = matrix(
    NA_real_,
    nrow = B,
    ncol = length(parameterNames),
    dimnames = list(NULL, parameterNames)
  )
  probabilityMatrix = matrix(
    NA_real_,
    nrow = B,
    ncol = length(probabilityNames),
    dimnames = list(NULL, probabilityNames)
  )
  failureMessages = rep(NA_character_, B)

  for (replicateIndex in seq_len(B)) {
    replicateFit = fitBayesianBootstrapReplicate(
      x = x,
      model = model,
      weights = weightDraws[replicateIndex, ],
      nterms = nterms
    )

    if (is.null(replicateFit$error)) {
      parameterMatrix[replicateIndex, ] = replicateFit$parameters
      probabilityMatrix[replicateIndex, ] = replicateFit$probabilities
    } else {
      failureMessages[replicateIndex] = replicateFit$error
    }
  }

  parameterReplicates = as.data.frame(parameterMatrix)
  probabilityReplicates = as.data.frame(probabilityMatrix)
  successful = complete.cases(parameterReplicates) &
    complete.cases(probabilityReplicates)
  nSuccessful = sum(successful)
  nFailed = B - nSuccessful

  diagnostics = list(
    nRequested = B,
    nSuccessful = nSuccessful,
    nFailed = nFailed,
    failureRate = nFailed / B,
    failures = data.frame(
      replicate = which(!successful),
      message = failureMessages[!successful],
      stringsAsFactors = FALSE
    )
  )

  parameterSummaries = summariseBayesianBootstrapReplicates(
    parameterReplicates,
    nameColumn = "parameter",
    level = level
  )
  probabilitySummaries = summariseBayesianBootstrapReplicates(
    probabilityReplicates,
    nameColumn = "term",
    level = level
  )

  newPsBayesianBootstrap(
    method = "rubin",
    parameters = parameterSummaries,
    probabilities = probabilitySummaries,
    replicates = list(
      parameters = parameterReplicates,
      probabilities = probabilityReplicates
    ),
    level = level,
    diagnostics = diagnostics,
    seed = seed,
    model = model$model
  )
}
