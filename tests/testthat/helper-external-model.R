externalZeroBasedSupport = function(n, type) {
  type = match.arg(type, c("P", "S"))
  if (identical(type, "P")) {
    return(n)
  }
  n - 1L
}

externalPoissonModel = function() {
  psModel(
    model = "poisson",
    parameterNames = "lambda",
    subclass = "externalPoissonModel",
    supportedEngines = c("numerical", "mcmc"),
    mleStart = c(lambda = 1),
    mleLower = c(lambda = sqrt(.Machine$double.eps))
  )
}

modelObservationData.externalPoissonModel = function(model, x, ...) {
  externalZeroBasedSupport(x$data$n, x$type)
}

modelProbabilities.externalPoissonModel = function(model,
                                                     parameters,
                                                     n,
                                                     type,
                                                     ...) {
  lambda = parameters[["lambda"]]
  support = externalZeroBasedSupport(n, type)
  values = vapply(
    support,
    function(value) {
      dpois(value, lambda = lambda)
    },
    numeric(length(lambda))
  )
  if (length(lambda) == 1L) {
    values = matrix(values, nrow = 1L)
  }
  colnames(values) = paste0(type, n)
  values
}

modelLogLikelihood.externalPoissonModel = function(model,
                                                    parameters,
                                                    data,
                                                    ...) {
  lambda = parameters[["lambda"]]
  observations = modelObservationData(model, data)
  sum(data$data$rn * dpois(observations, lambda = lambda, log = TRUE))
}

registerExternalPoissonMethods = function() {
  namespace = asNamespace("fitPS")
  registerS3method(
    "modelObservationData",
    "externalPoissonModel",
    modelObservationData.externalPoissonModel,
    envir = namespace
  )
  registerS3method(
    "modelProbabilities",
    "externalPoissonModel",
    modelProbabilities.externalPoissonModel,
    envir = namespace
  )
  registerS3method(
    "modelLogLikelihood",
    "externalPoissonModel",
    modelLogLikelihood.externalPoissonModel,
    envir = namespace
  )
}

registerExternalPoissonMethods()


modelLogPrior.externalPoissonModel = function(model, parameters, prior, ...) {
  lambda = parameters[["lambda"]]
  if (!is.list(prior) ||
      !all(c("shape", "rate") %in% names(prior)) ||
      any(!is.finite(c(prior$shape, prior$rate))) ||
      prior$shape <= 0 || prior$rate <= 0) {
    stop("Poisson prior must contain positive finite shape and rate values")
  }
  dgamma(lambda, shape = prior$shape, rate = prior$rate, log = TRUE)
}

modelBayesControl.externalPoissonModel = function(model, x, engine, prior, ...) {
  observations = modelObservationData(model, x)
  start = weighted.mean(observations, x$data$rn)
  if (!is.finite(start) || start <= 0) {
    start = 1
  }
  list(
    start = c(lambda = start),
    lower = c(lambda = 0),
    upper = c(lambda = Inf)
  )
}

modelToUnconstrained.externalPoissonModel = function(model, parameters, ...) {
  parameters = validateExternalPoissonParameters(model, parameters)
  c(lambda = log(parameters[["lambda"]]))
}

modelFromUnconstrained.externalPoissonModel = function(model, unconstrained, ...) {
  unconstrained = validateExternalPoissonParameters(model, unconstrained, positive = FALSE)
  c(lambda = exp(unconstrained[["lambda"]]))
}

modelLogJacobian.externalPoissonModel = function(model, unconstrained, ...) {
  unconstrained = validateExternalPoissonParameters(model, unconstrained, positive = FALSE)
  unname(unconstrained[["lambda"]])
}

validateExternalPoissonParameters = function(model, value, positive = TRUE) {
  if (!is.numeric(value) || length(value) != 1L ||
      is.null(names(value)) || !identical(names(value), "lambda") ||
      !is.finite(value[["lambda"]])) {
    stop("value must contain one finite named lambda parameter")
  }
  if (positive && value[["lambda"]] <= 0) {
    stop("lambda must be positive")
  }
  value
}

registerExternalPoissonBayesMethods = function() {
  namespace = asNamespace("fitPS")
  methods = c(
    modelLogPrior = "modelLogPrior.externalPoissonModel",
    modelBayesControl = "modelBayesControl.externalPoissonModel",
    modelToUnconstrained = "modelToUnconstrained.externalPoissonModel",
    modelFromUnconstrained = "modelFromUnconstrained.externalPoissonModel",
    modelLogJacobian = "modelLogJacobian.externalPoissonModel"
  )
  for (generic in names(methods)) {
    registerS3method(
      generic,
      "externalPoissonModel",
      get(methods[[generic]], mode = "function"),
      envir = namespace
    )
  }
}

registerExternalPoissonBayesMethods()
