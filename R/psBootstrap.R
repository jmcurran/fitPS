#' Construct a fitPS bootstrap object
#'
#' @param method Bootstrap method description.
#' @param parameters Data frame of bootstrap parameter summaries.
#' @param probabilities Data frame of bootstrap probability summaries.
#' @param replicates Data frame containing the raw bootstrap parameter
#'   replicates. Failed fits are retained as rows containing missing values.
#' @param level Confidence interval level used for the summaries.
#' @param diagnostics Bootstrap diagnostics.
#'
#' @return An object of class `psBootstrap`.
#'
#' @importFrom doParallel registerDoParallel
#' @import foreach
#' @importFrom iterators iter
#' @importFrom parallel detectCores makeCluster parLapply stopCluster
#' @importFrom pbapply pbapply pblapply pboptions
#' @importFrom stats complete.cases quantile sd
#' @keywords internal
newPsBootstrap = function(method,
                           parameters,
                           probabilities,
                           replicates,
                           level = 0.95,
                           diagnostics = NULL) {
  if (!is.character(method) || length(method) != 1L || !nzchar(method)) {
    stop("method must be one non-empty character value")
  }

  requiredParameterColumns = c(
    "parameter", "estimate", "sd", "lower", "upper", "level"
  )
  if (!is.data.frame(parameters) ||
      !all(requiredParameterColumns %in% names(parameters))) {
    stop("parameters must contain bootstrap parameter summaries")
  }

  requiredProbabilityColumns = c(
    "term", "estimate", "sd", "lower", "upper", "level",
    "bootstrapMethod"
  )
  if (!is.data.frame(probabilities) ||
      !all(requiredProbabilityColumns %in% names(probabilities))) {
    stop("probabilities must contain bootstrap probability summaries")
  }

  if (!is.data.frame(replicates)) {
    stop("replicates must be a data frame")
  }

  result = list(
    method = method,
    parameters = parameters,
    probabilities = probabilities,
    replicates = replicates,
    level = level,
    diagnostics = diagnostics
  )
  class(result) = "psBootstrap"
  result
}

#' Validate a bootstrap confidence level.
#'
#' @param level Probability level for intervals or summaries.
#' @return The validated level, invisibly.
#' @keywords internal
#' @noRd
validateBootstrapLevel = function(level) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("level must be one finite number strictly between 0 and 1")
  }
  invisible(level)
}

#' Compute equal-tail bootstrap quantiles at a requested confidence level.
#'
#' @param x An input object or numeric vector required by the helper.
#' @param level Probability level for intervals or summaries.
#' @return Numeric lower and upper quantiles.
#' @keywords internal
#' @noRd
bootstrapQuantiles = function(x, level) {
  alpha = (1 - level) / 2
  unname(quantile(
    x,
    probs = c(alpha, 1 - alpha),
    names = FALSE,
    type = 7
  ))
}

#' Summarise bootstrap parameter replicates with estimates, uncertainty, and intervals.
#'
#' @param replicates Bootstrap parameter replicate matrix or data frame.
#' @param level Probability level for intervals or summaries.
#' @return A data frame of bootstrap parameter summaries.
#' @keywords internal
#' @noRd
summariseBootstrapParameters = function(replicates, level) {
  validateBootstrapLevel(level)

  summaries = lapply(names(replicates), function(parameterName) {
    values = replicates[[parameterName]]
    values = values[is.finite(values)]
    bounds = bootstrapQuantiles(values, level = level)

    data.frame(
      parameter = parameterName,
      estimate = mean(values),
      sd = sd(values),
      lower = bounds[1L],
      upper = bounds[2L],
      level = level,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, summaries)
}

#' Summarise bootstrap fitted-probability replicates.
#'
#' @param probabilities Bootstrap or posterior fitted-probability values.
#' @param level Probability level for intervals or summaries.
#' @param method Character string identifying a fitting or posterior method.
#' @return A data frame of bootstrap probability summaries.
#' @keywords internal
#' @noRd
summariseBootstrapProbabilities = function(probabilities,
                                            level,
                                            method = "nonparametric") {
  validateBootstrapLevel(level)

  summaries = lapply(seq_len(ncol(probabilities)), function(columnIndex) {
    values = probabilities[, columnIndex]
    values = values[is.finite(values)]
    bounds = bootstrapQuantiles(values, level = level)

    data.frame(
      term = colnames(probabilities)[columnIndex],
      estimate = mean(values),
      sd = sd(values),
      lower = bounds[1L],
      upper = bounds[2L],
      level = level,
      bootstrapMethod = method,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, summaries)
}

#' Fit the requested model to one bootstrap resample.
#'
#' @param x An input object or numeric vector required by the helper.
#' @param model Model identifier, typically `"zeta"` or `"ziz"`.
#' @return A fitted zeta or ZIZ model.
#' @keywords internal
#' @noRd
fitBootstrapSample = function(x, model) {
  tryCatch(
    {
      if (model == "zeta") {
        fit = fitDist(x)
        return(c(shape = unname(fit$shape)))
      }

      fit = fitZIDist(x)
      c(pi = unname(fit$pi), shape = unname(fit$shape))
    },
    error = function(e) {
      if (model == "zeta") {
        return(c(shape = NA_real_))
      }
      c(pi = NA_real_, shape = NA_real_)
    }
  )
}

#' Generate model-parameter bootstrap replicates with optional parallel execution.
#'
#' @param x An input object or numeric vector required by the helper.
#' @param B Number of bootstrap replicates.
#' @param model Model identifier, typically `"zeta"` or `"ziz"`.
#' @param seed Optional random-number seed.
#' @param silent Logical; suppress progress output when `TRUE`.
#' @param parallel Logical; use parallel computation when supported.
#' @param progressBar Logical; display progress information when supported.
#' @param pbopts Progress-bar options passed to `pbapply` helpers.
#' @return A collection of bootstrap parameter replicates.
#' @keywords internal
#' @noRd
bootstrapParameterReplicates = function(x,
                                         B = 2000,
                                         model = c("zeta", "ziz"),
                                         seed = NULL,
                                         silent = FALSE,
                                         parallel = TRUE,
                                         progressBar = FALSE,
                                         pbopts = list(type = "txt")) {
  if (!inherits(x, "psData")) {
    stop("x must be an object of class psData")
  }
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) ||
      B < 1 || B != floor(B)) {
    stop("B must be one positive integer")
  }
  B = as.integer(B)
  model = match.arg(model)

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("seed must be NULL or one finite numeric value")
    }
    set.seed(as.integer(seed))
  }

  yvals = rep(x$data$n, x$data$rn)
  sampleSize = length(yvals)
  bootstrapMatrix = matrix(
    sample(yvals, sampleSize * B, replace = TRUE),
    nrow = B
  )

  toPsData = function(y, type) {
    tbl = table(y)
    result = list(
      data = data.frame(
        n = as.numeric(names(tbl)),
        rn = as.vector(tbl)
      ),
      type = type
    )
    class(result) = "psData"
    result
  }

  if (!silent) {
    cat("Creating bootstrapped data sets\n")
  }

  if (progressBar) {
    oldPbOptions = do.call(
      get("pboptions", asNamespace("pbapply")),
      pbopts
    )
    on.exit(pboptions(oldPbOptions), add = TRUE)
  }

  bootstrapData = if (parallel) {
    nCores = detectCores()
    cluster = makeCluster(nCores, setup_strategy = "sequential")
    on.exit(stopCluster(cluster), add = TRUE)
    registerDoParallel(cluster)

    if (progressBar) {
      pbapply(
        X = bootstrapMatrix,
        MARGIN = 1,
        FUN = toPsData,
        type = x$type,
        cl = cluster
      )
    } else {
      foreach(row = iter(bootstrapMatrix, by = "row")) %dopar% {
        toPsData(row, type = x$type)
      }
    }
  } else if (progressBar) {
    pbapply(
      X = bootstrapMatrix,
      MARGIN = 1,
      FUN = toPsData,
      type = x$type
    )
  } else {
    apply(bootstrapMatrix, 1, toPsData, type = x$type)
  }

  if (!silent) {
    cat("Estimating parameters for each bootstrapped data set\n")
  }

  fittedValues = if (parallel) {
    if (progressBar) {
      pblapply(
        X = bootstrapData,
        FUN = fitBootstrapSample,
        model = model,
        cl = cluster
      )
    } else {
      parLapply(
        cl = cluster,
        X = bootstrapData,
        fun = fitBootstrapSample,
        model = model
      )
    }
  } else if (progressBar) {
    pblapply(
      X = bootstrapData,
      FUN = fitBootstrapSample,
      model = model
    )
  } else {
    lapply(
      bootstrapData,
      fitBootstrapSample,
      model = model
    )
  }

  replicates = as.data.frame(do.call(rbind, fittedValues))
  expectedNames = if (model == "zeta") "shape" else c("pi", "shape")
  names(replicates) = expectedNames

  successful = complete.cases(replicates)
  nSuccessful = sum(successful)
  nFailed = B - nSuccessful

  if (nSuccessful == 0L) {
    stop("All bootstrap model fits failed")
  }

  diagnostics = list(
    B = B,
    nSuccessful = nSuccessful,
    nFailed = nFailed,
    failureRate = nFailed / B,
    seed = seed
  )

  list(
    replicates = replicates,
    successful = successful,
    diagnostics = diagnostics
  )
}

#' Construct a bootstrap summary for an existing fitted P/S model.
#'
#' @param object Posterior or fitted object.
#' @param B Number of bootstrap replicates.
#' @param level Probability level for intervals or summaries.
#' @param seed Optional random-number seed.
#' @param silent Logical; suppress progress output when `TRUE`.
#' @param parallel Logical; use parallel computation when supported.
#' @param progressBar Logical; display progress information when supported.
#' @param pbopts Progress-bar options passed to `pbapply` helpers.
#' @return A `psBootstrap` object.
#' @keywords internal
#' @noRd
bootstrapPsFit = function(object,
                           B = 2000,
                           level = 0.95,
                           seed = NULL,
                           silent = FALSE,
                           parallel = TRUE,
                           progressBar = FALSE,
                           pbopts = list(type = "txt")) {
  if (!inherits(object, "psFit")) {
    stop("object must be an object of class psFit")
  }
  if (!identical(object$method, "mle")) {
    stop("bootstrap inference currently requires an MLE psFit object")
  }
  validateBootstrapLevel(level)

  bootstrapResult = bootstrapParameterReplicates(
    x = object$psData,
    B = B,
    model = object$model,
    seed = seed,
    silent = silent,
    parallel = parallel,
    progressBar = progressBar,
    pbopts = pbopts
  )

  successfulReplicates = bootstrapResult$replicates[
    bootstrapResult$successful,
    ,
    drop = FALSE
  ]

  parameterSummary = summariseBootstrapParameters(
    replicates = successfulReplicates,
    level = level
  )

  nTerms = length(object$fitted)
  termIndices = if (object$psData$type == "P") {
    seq.int(0L, nTerms - 1L)
  } else {
    seq_len(nTerms)
  }

  probabilityReplicates = if (object$model == "ziz") {
    zizProbabilities(
      pi = successfulReplicates$pi,
      shape = successfulReplicates$shape,
      n = termIndices,
      type = object$psData$type
    )
  } else {
    zetaProbabilities(
      shape = successfulReplicates$shape,
      n = termIndices,
      type = object$psData$type
    )
  }

  probabilitySummary = summariseBootstrapProbabilities(
    probabilities = probabilityReplicates,
    level = level,
    method = "nonparametric"
  )

  object$bootstrap = newPsBootstrap(
    method = "nonparametric",
    parameters = parameterSummary,
    probabilities = probabilitySummary,
    replicates = bootstrapResult$replicates,
    level = level,
    diagnostics = bootstrapResult$diagnostics
  )
  object
}

#' Print a fitPS bootstrap object
#'
#' @param x An object of class `psBootstrap`.
#' @param nterms Optional number of leading probability summaries to print.
#'   By default, at most the first 10 are shown.
#' @param ... Additional arguments passed to `print.data.frame()`.
#'
#' @return The bootstrap object, invisibly.
#' @export
print.psBootstrap = function(x, nterms = NULL, ...) {
  if (is.null(nterms)) {
    nterms = min(10L, nrow(x$probabilities))
  }
  probabilities = bootstrapProbs(x, n = nterms)

  cat("fitPS bootstrap distribution\n")
  cat("Method:", x$method, "\n\n")
  cat("Parameter summaries:\n")
  print(x$parameters, row.names = FALSE, ...)
  cat("\nBootstrap probability summaries:\n")
  print(probabilities, row.names = FALSE, ...)
  if (nrow(probabilities) < nrow(x$probabilities)) {
    cat(
      "\nShowing", nrow(probabilities), "of", nrow(x$probabilities),
      "stored probability summaries.\n"
    )
  }

  if (!is.null(x$diagnostics)) {
    cat("\nDiagnostics:\n")
    print(x$diagnostics)
  }

  invisible(x)
}

#' Summarise a fitPS bootstrap object
#'
#' @param object An object of class `psBootstrap`.
#' @param nterms Optional number of leading probability summaries to include.
#'   If `NULL`, all stored summaries are included.
#' @param ... Additional arguments retained for S3 compatibility.
#'
#' @return An object of class `summary.psBootstrap`.
#' @export
summary.psBootstrap = function(object, nterms = NULL, ...) {
  probabilities = if (is.null(nterms)) {
    object$probabilities
  } else {
    bootstrapProbs(object, n = nterms)
  }

  result = list(
    method = object$method,
    parameters = object$parameters,
    probabilities = probabilities,
    level = object$level,
    diagnostics = object$diagnostics
  )
  class(result) = "summary.psBootstrap"
  result
}

#' @describeIn summary.psBootstrap Print a summarized fitPS bootstrap object.
#' @param x An object of class `summary.psBootstrap`.
#' @export
print.summary.psBootstrap = function(x, ...) {
  cat("Summary of fitPS bootstrap distribution\n")
  cat("Method:", x$method, "\n\n")
  cat("Parameter summaries:\n")
  print(x$parameters, row.names = FALSE, ...)
  cat("\nBootstrap probability summaries:\n")
  print(x$probabilities, row.names = FALSE, ...)

  if (!is.null(x$diagnostics)) {
    cat("\nDiagnostics:\n")
    print(x$diagnostics)
  }

  invisible(x)
}

