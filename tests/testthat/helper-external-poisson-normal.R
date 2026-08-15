externalPoissonNormalModel = function() {
  psModel(
    model = "poissonNormal",
    parameterNames = c("mu", "sigma"),
    subclass = "externalPoissonNormalModel"
  )
}

modelObservationData.externalPoissonNormalModel = function(model, x, ...) {
  x$data$n
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
  mu = parameters[["mu"]]
  sigma = parameters[["sigma"]]
  values = externalPoissonNormalProbability(n, mu = mu, sigma = sigma)
  values = matrix(values, nrow = 1L)
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
