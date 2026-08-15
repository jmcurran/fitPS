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
  muStart = weighted.mean(observations, weights)
  varianceStart = weighted.mean((observations - muStart)^2, weights)
  lowerVariance = muStart * 1.05
  upperVariance = muStart * 1.95
  varianceStart = min(max(varianceStart, lowerVariance), upperVariance)

  list(
    start = c(mu = muStart, sigma = sqrt(varianceStart)),
    lower = c(mu = sqrt(.Machine$double.eps), sigma = sqrt(.Machine$double.eps)),
    upper = c(mu = Inf, sigma = Inf)
  )
}

externalPoissonNormalComponentRates = function(mu, sigma) {
  variance = sigma^2
  lambdaTwo = (variance - mu) / 2
  lambdaOne = 2 * mu - variance

  if (!is.finite(mu) || !is.finite(sigma) || mu <= 0 || sigma <= 0 ||
      lambdaOne < 0 || lambdaTwo < 0) {
    return(NULL)
  }

  c(lambdaOne = lambdaOne, lambdaTwo = lambdaTwo)
}

externalPoissonNormalProbability = function(n, mu, sigma) {
  rates = externalPoissonNormalComponentRates(mu, sigma)
  if (is.null(rates)) {
    return(rep(NaN, length(n)))
  }

  vapply(n, function(value) {
    if (!is.finite(value) || value < 0 || value != floor(value)) {
      return(0)
    }

    doubleCount = 0:floor(value / 2)
    singleCount = value - 2 * doubleCount
    sum(
      dpois(singleCount, lambda = rates[["lambdaOne"]]) *
        dpois(doubleCount, lambda = rates[["lambdaTwo"]])
    )
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
