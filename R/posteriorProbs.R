#' Extract posterior probability summaries
#'
#' Extract posterior summaries for fitted P or S probabilities from a Bayesian
#' fitPS model or from its associated posterior object.
#'
#' @param object A Bayesian `psFit` object or a `psPosterior` object.
#' @param n `NULL`, a single number of leading probabilities to return, or a
#'   vector of P or S indices to return.
#' @param ... Additional arguments passed to methods.
#'
#' @return A data frame containing posterior probability summaries.
#'
#' @export
posteriorProbs = function(object, ...) {
  UseMethod("posteriorProbs")
}

#' @rdname posteriorProbs
#' @export
posteriorProbs.psFit = function(object, n = NULL, ...) {
  if (!identical(object$method, "bayes") ||
      !inherits(object$posterior, "psPosterior")) {
    stop("posteriorProbs() is only available for Bayesian psFit objects")
  }

  posteriorProbs(object$posterior, n = n, ...)
}

#' @rdname posteriorProbs
#' @export
posteriorProbs.psPosterior = function(object, n = NULL, ...) {
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
        "posterior summaries are not available for: ",
        paste(missingTerms, collapse = ", ")
      )
    }
  }

  probabilities[selectedRows, , drop = FALSE]
}
