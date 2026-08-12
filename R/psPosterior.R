#' Construct a fitPS posterior object
#'
#' @param method Posterior approximation method.
#' @param parameters Data frame of posterior parameter summaries.
#' @param probabilities Data frame of posterior probability summaries.
#' @param representation Engine-specific posterior representation.
#' @param level Credible interval level used for probability summaries.
#' @param diagnostics Optional posterior diagnostics.
#' @param model Model identifier associated with the posterior.
#'
#' @return An object of class `psPosterior`.
#'
#' @keywords internal
newPsPosterior = function(method,
                           parameters,
                           probabilities,
                           representation,
                           level = 0.95,
                           diagnostics = NULL,
                           model = NULL) {
  if (!is.character(method) || length(method) != 1L || !nzchar(method)) {
    stop("method must be one non-empty character value")
  }

  if (!is.data.frame(parameters) ||
      !all(c("parameter", "estimate", "sd") %in% names(parameters))) {
    stop("parameters must contain parameter, estimate, and sd columns")
  }

  if (!is.data.frame(probabilities) ||
      !all(c("term", "estimate", "sd", "lower", "upper", "level") %in%
        names(probabilities))) {
    stop("probabilities must contain posterior probability summaries")
  }

  result = list(
    method = method,
    parameters = parameters,
    probabilities = probabilities,
    representation = representation,
    level = level,
    diagnostics = diagnostics,
    model = model
  )
  class(result) = "psPosterior"
  result
}

makeZizParameterSummary = function(par, varCov) {
  par = unname(par)
  names(par) = c("pi", "shape")
  standardErrors = sqrt(pmax(0, diag(varCov)))

  data.frame(
    parameter = names(par),
    estimate = unname(par),
    sd = unname(standardErrors),
    stringsAsFactors = FALSE
  )
}

attachZizPosterior = function(result,
                               probabilities,
                               representation,
                               diagnostics = NULL,
                               level = 0.95) {
  posterior = newPsPosterior(
    method = result$posteriorMethod,
    parameters = makeZizParameterSummary(
      par = c(pi = result$pi, shape = result$shape),
      varCov = result$var.cov
    ),
    probabilities = probabilities,
    representation = representation,
    level = level,
    diagnostics = diagnostics,
    model = "ziz"
  )

  result$posterior = posterior
  result$posteriorProbs = posterior$probabilities
  result
}

#' Print a fitPS posterior object
#'
#' @param x An object of class `psPosterior`.
#' @param nterms Optional number of leading probability summaries to print.
#'   By default, at most the first 10 are shown.
#' @param ... Additional arguments passed to `print.data.frame()`.
#'
#' @return The posterior object, invisibly.
#' @export
print.psPosterior = function(x, nterms = NULL, ...) {
  if (is.null(nterms)) {
    nterms = min(10L, nrow(x$probabilities))
  }
  probabilities = posteriorProbs(x, n = nterms)

  cat("fitPS posterior approximation\n")
  cat("Method:", x$method, "\n\n")
  cat("Parameter summaries:\n")
  print(x$parameters, row.names = FALSE, ...)
  cat("\nPosterior probability summaries:\n")
  print(probabilities, row.names = FALSE, ...)
  if (nrow(probabilities) < nrow(x$probabilities)) {
    cat(
      "\nShowing", nrow(probabilities), "of", nrow(x$probabilities),
      "stored probability summaries.\n"
    )
  }
  invisible(x)
}

#' Summarise a fitPS posterior object
#'
#' @param object An object of class `psPosterior`.
#' @param nterms Optional number of leading probability summaries to include.
#'   If `NULL`, all stored summaries are included.
#' @param inflationEpsilon Practical threshold for the inflation parameter.
#'   For zero-inflated posteriors, the summary reports
#'   `Pr(pi < inflationEpsilon | data)`. The default is 0.01.
#' @param ... Additional arguments retained for S3 compatibility.
#'
#' @return An object of class `summary.psPosterior`.
#' @export
summary.psPosterior = function(object,
                                nterms = NULL,
                                inflationEpsilon = 0.01,
                                ...) {
  probabilities = if (is.null(nterms)) {
    object$probabilities
  } else {
    posteriorProbs(object, n = nterms)
  }

  result = list(
    method = object$method,
    parameters = object$parameters,
    probabilities = probabilities,
    level = object$level,
    diagnostics = object$diagnostics,
    inflation = if (identical(object$model, "ziz")) {
      posteriorInflation(object, epsilon = inflationEpsilon)
    } else {
      NULL
    }
  )
  class(result) = "summary.psPosterior"
  result
}

#' @export
print.summary.psPosterior = function(x, ...) {
  cat("Summary of fitPS posterior approximation\n")
  cat("Method:", x$method, "\n\n")
  cat("Parameter summaries:\n")
  print(x$parameters, row.names = FALSE, ...)
  cat("\nPosterior probability summaries:\n")
  print(x$probabilities, row.names = FALSE, ...)

  if (!is.null(x$inflation)) {
    cat("\nInflation diagnostic:\n")
    cat(
      "Pr(pi < ", format(x$inflation$epsilon), " | data) = ",
      format(x$inflation$probBelow, digits = 4), "\n",
      sep = ""
    )
  }

  if (!is.null(x$diagnostics)) {
    cat("\nDiagnostics:\n")
    print(x$diagnostics)
  }

  invisible(x)
}

#' Extract posterior mean probabilities
#'
#' @param object An object of class `psPosterior`.
#' @param ... Additional arguments retained for S3 compatibility.
#'
#' @return A named numeric vector of posterior mean probabilities.
#' @export
fitted.psPosterior = function(object, ...) {
  estimates = object$probabilities$estimate
  names(estimates) = object$probabilities$term
  estimates
}
