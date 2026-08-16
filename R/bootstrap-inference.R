#' Compute a bootstrap distribution for a fitted fitPS model
#'
#' `bootstrapFit()` is a deprecated compatibility wrapper. New code should use
#' [fit()] with `method = "bootstrap"`. The wrapper computes and attaches a
#' nonparametric bootstrap distribution to a maximum-likelihood `psFit` object.
#'
#' @param object A maximum likelihood \code{psFit} object.
#' @param B Number of bootstrap replicates.
#' @param level Confidence level used for percentile intervals.
#' @param seed Optional random seed for reproducible resampling.
#' @param silent Logical; suppress progress messages when \code{TRUE}.
#' @param parallel Logical; use parallel fitting when \code{TRUE}.
#' @param progressBar Logical; display a progress bar when \code{TRUE}.
#' @param pbopts Options passed to \code{pbapply::pboptions()} when progress
#'   bars are requested.
#' @param ... Additional arguments passed to methods.
#'
#' @return The input \code{psFit} object with a \code{psBootstrap} object
#'   attached as \code{bootstrap}.
#'
#' @details
#' The bootstrap is an observation-level nonparametric bootstrap. Each
#' successful parameter replicate is transformed into the corresponding P or S
#' probabilities before probability summaries are calculated. Consequently,
#' the reported bootstrap mean probabilities are averages of transformed
#' bootstrap replicates, not probabilities evaluated at average bootstrap
#' parameter values.
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
#' bootstrapProbs(fit)
#' }
#' @export
bootstrapFit = function(object, ...) {
  signalDeprecatedFitter(
    "bootstrapFit",
    "fit(x, model = ..., method = \"bootstrap\")"
  )
  UseMethod("bootstrapFit")
}

#' @rdname bootstrapFit
#' @export
bootstrapFit.psFit = function(object,
                               B = 2000,
                               level = 0.95,
                               seed = NULL,
                               silent = FALSE,
                               parallel = TRUE,
                               progressBar = FALSE,
                               pbopts = list(type = "txt"),
                               ...) {
  bootstrapPsFit(
    object = object,
    B = B,
    level = level,
    seed = seed,
    silent = silent,
    parallel = parallel,
    progressBar = progressBar,
    pbopts = pbopts
  )
}
#' Extract bootstrap probability summaries
#'
#' Extract bootstrap summaries for fitted P or S probabilities from a
#' frequentist fitPS model with an attached bootstrap distribution or directly
#' from its \code{psBootstrap} object.
#'
#' @param object A \code{psFit} object with an attached \code{psBootstrap}
#'   object, or a \code{psBootstrap} object.
#' @param n \code{NULL}, a single number of leading probabilities to return,
#'   or a vector of P or S indices to return.
#' @param ... Additional arguments passed to methods.
#'
#' @return A data frame containing bootstrap probability summaries.
#'
#' @details
#' The \code{estimate} column is the bootstrap mean of the probability across
#' successful bootstrap replicates. It is not a Bayesian posterior mean and is
#' generally not identical to the plug-in probability evaluated at the MLE.
#' The reported bounds are percentile bootstrap confidence intervals.
#'
#' @export
bootstrapProbs = function(object, ...) {
  UseMethod("bootstrapProbs")
}

#' @rdname bootstrapProbs
#' @export
bootstrapProbs.psFit = function(object, n = NULL, ...) {
  if (!inherits(object$bootstrap, "psBootstrap")) {
    stop(
      "bootstrapProbs() requires a psFit object with an attached bootstrap; ",
      "use fit(..., method = 'bootstrap')"
    )
  }

  bootstrapProbs(object$bootstrap, n = n, ...)
}

#' @rdname bootstrapProbs
#' @export
bootstrapProbs.psBootstrap = function(object, n = NULL, ...) {
  probabilities = object$probabilities

  if (is.null(n)) {
    return(probabilities)
  }

  if (!is.numeric(n) || length(n) < 1L || any(!is.finite(n))) {
    stop("n must be NULL or a finite numeric vector")
  }

  if (any(abs(n - floor(n)) > .Machine$double.eps^0.5)) {
    stop("n must contain integer values")
  }

  n = as.integer(n)
  termPrefix = substring(probabilities$term[1L], 1L, 1L)
  minimumIndex = if (identical(termPrefix, "P")) 0L else 1L

  if (length(n) == 1L) {
    if (n < 1L) {
      stop("a single n value must be at least 1")
    }

    selectedRows = seq_len(min(n, nrow(probabilities)))
  } else {
    if (any(n < minimumIndex)) {
      stop(
        if (identical(termPrefix, "P")) {
          "P indices must be at least 0"
        } else {
          "S indices must be at least 1"
        }
      )
    }

    requestedTerms = paste0(termPrefix, n)
    selectedRows = match(requestedTerms, probabilities$term)

    if (anyNA(selectedRows)) {
      missingTerms = requestedTerms[is.na(selectedRows)]
      stop(
        "bootstrap summaries are not available for: ",
        paste(missingTerms, collapse = ", ")
      )
    }
  }

  probabilities[selectedRows, , drop = FALSE]
}
