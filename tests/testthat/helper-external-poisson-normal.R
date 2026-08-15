externalPoissonNormalModel = function() {
  psModel(
    model = "poissonNormal",
    parameterNames = c("mu", "sigma"),
    subclass = "externalPoissonNormalModel",
    supportedEngines = "mcmc"
  )
}

modelObservationData.externalPoissonNormalModel = function(model, x, ...) {
  externalZeroBasedSupport(x$data$n, x$type)
}

modelMleControl.externalPoissonNormalModel = function(model, x, ...) {
  observations = modelObservationData(model, x)
  weights = x$data$rn
  meanStart = weighted.mean(observations, weights)
  varianceStart = weighted.mean((observations - meanStart)^2, weights)

  if (!is.finite(meanStart) || meanStart <= 0) {
    meanStart = 1
  }
  if (!is.finite(varianceStart)) {
    varianceStart = meanStart
  }

  extraVariance = max(
    varianceStart - meanStart,
    meanStart^2 * 1e-6
  )
  sigmaSquaredStart = log1p(extraVariance / meanStart^2)
  sigmaStart = sqrt(max(sigmaSquaredStart, 1e-6))
  muStart = log(meanStart) - sigmaSquaredStart / 2

  list(
    start = c(mu = muStart, sigma = sigmaStart),
    lower = c(mu = -20, sigma = sqrt(.Machine$double.eps)),
    upper = c(mu = 20, sigma = 5)
  )
}

externalPoissonNormalProbability = function(n, mu, sigma) {
  if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0) {
    return(rep(NaN, length(n)))
  }

  vapply(n, function(value) {
    if (!is.finite(value) || value < 0 || value != floor(value)) {
      return(0)
    }

    integrand = function(z) {
      dpois(value, lambda = exp(z)) * dnorm(z, mean = mu, sd = sigma)
    }

    result = integrate(
      integrand,
      lower = -Inf,
      upper = Inf,
      rel.tol = 1e-8,
      subdivisions = 200L,
      stop.on.error = FALSE
    )

    if (!identical(result$message, "OK")) {
      return(NaN)
    }

    result$value
  }, numeric(1L))
}

modelProbabilities.externalPoissonNormalModel = function(model,
                                                          parameters,
                                                          n,
                                                          type,
                                                          ...) {
  parameterFrame = as.data.frame(parameters)
  if (!all(c("mu", "sigma") %in% names(parameterFrame))) {
    stop("parameters must contain mu and sigma")
  }
  support = externalZeroBasedSupport(n, type)
  values = vapply(seq_len(nrow(parameterFrame)), function(row) {
    externalPoissonNormalProbability(
      support,
      mu = parameterFrame$mu[row],
      sigma = parameterFrame$sigma[row]
    )
  }, numeric(length(support)))
  values = t(values)
  colnames(values) = paste0(type, n)
  values
}

modelLogLikelihood.externalPoissonNormalModel = function(model,
                                                           parameters,
                                                           data,
                                                           ...) {
  mu = parameters[["mu"]]
  sigma = parameters[["sigma"]]
  observations = modelObservationData(model, data)
  probabilities = externalPoissonNormalProbability(
    observations,
    mu = mu,
    sigma = sigma
  )

  if (any(!is.finite(probabilities)) || any(probabilities <= 0)) {
    return(-Inf)
  }

  sum(data$data$rn * log(probabilities))
}

registerExternalPoissonNormalMethods = function() {
  namespace = asNamespace("fitPS")
  registerS3method(
    "modelObservationData",
    "externalPoissonNormalModel",
    modelObservationData.externalPoissonNormalModel,
    envir = namespace
  )
  registerS3method(
    "modelMleControl",
    "externalPoissonNormalModel",
    modelMleControl.externalPoissonNormalModel,
    envir = namespace
  )
  registerS3method(
    "modelProbabilities",
    "externalPoissonNormalModel",
    modelProbabilities.externalPoissonNormalModel,
    envir = namespace
  )
  registerS3method(
    "modelLogLikelihood",
    "externalPoissonNormalModel",
    modelLogLikelihood.externalPoissonNormalModel,
    envir = namespace
  )
}

registerExternalPoissonNormalMethods()


modelLogPrior.externalPoissonNormalModel = function(model, parameters, prior, ...) {
  mu = parameters[["mu"]]
  sigma = parameters[["sigma"]]
  required = c("muMean", "muSd", "sigmaScale")
  if (!is.list(prior) || !all(required %in% names(prior)) ||
      any(!is.finite(unlist(prior[required], use.names = FALSE))) ||
      prior$muSd <= 0 || prior$sigmaScale <= 0) {
    stop("Poisson-normal prior must contain muMean, positive muSd, and positive sigmaScale")
  }
  if (!is.finite(sigma) || sigma <= 0) {
    return(-Inf)
  }
  dnorm(mu, mean = prior$muMean, sd = prior$muSd, log = TRUE) +
    log(2) + dnorm(sigma, mean = 0, sd = prior$sigmaScale, log = TRUE)
}

modelBayesControl.externalPoissonNormalModel = function(model, x, engine, prior, ...) {
  control = modelMleControl(model, x)
  list(start = control$start)
}

modelToWorking.externalPoissonNormalModel = function(model, parameters, ...) {
  parameters = validateExternalPoissonNormalParameters(model, parameters)
  c(mu = parameters[["mu"]], sigma = log(parameters[["sigma"]]))
}

modelFromWorking.externalPoissonNormalModel = function(model, working, ...) {
  working = validateExternalPoissonNormalParameters(model, working, workingScale = TRUE)
  c(mu = working[["mu"]], sigma = exp(working[["sigma"]]))
}

modelWorkingLogJacobian.externalPoissonNormalModel = function(model, working, ...) {
  working = validateExternalPoissonNormalParameters(model, working, workingScale = TRUE)
  unname(working[["sigma"]])
}

validateExternalPoissonNormalParameters = function(model, value, workingScale = FALSE) {
  if (!is.numeric(value) || length(value) != 2L ||
      is.null(names(value)) || !setequal(names(value), c("mu", "sigma"))) {
    stop("value must be a named numeric vector containing mu and sigma")
  }
  value = value[c("mu", "sigma")]
  if (any(!is.finite(value))) {
    stop("mu and sigma must be finite")
  }
  if (!workingScale && value[["sigma"]] <= 0) {
    stop("sigma must be positive")
  }
  value
}

registerExternalPoissonNormalBayesMethods = function() {
  namespace = asNamespace("fitPS")
  methods = c(
    modelLogPrior = "modelLogPrior.externalPoissonNormalModel",
    modelBayesControl = "modelBayesControl.externalPoissonNormalModel",
    modelToWorking = "modelToWorking.externalPoissonNormalModel",
    modelFromWorking = "modelFromWorking.externalPoissonNormalModel",
    modelWorkingLogJacobian = "modelWorkingLogJacobian.externalPoissonNormalModel"
  )
  for (generic in names(methods)) {
    registerS3method(
      generic,
      "externalPoissonNormalModel",
      get(methods[[generic]], mode = "function"),
      envir = namespace
    )
  }
}

registerExternalPoissonNormalBayesMethods()
