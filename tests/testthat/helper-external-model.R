externalPoissonModel = function() {
  psModel(
    model = "poisson",
    parameterNames = "lambda",
    subclass = "externalPoissonModel",
    mleStart = c(lambda = 1),
    mleLower = c(lambda = sqrt(.Machine$double.eps))
  )
}

modelObservationData.externalPoissonModel = function(model, x, ...) {
  x$data$n
}

modelProbabilities.externalPoissonModel = function(model,
                                                     parameters,
                                                     n,
                                                     type,
                                                     ...) {
  lambda = parameters[["lambda"]]
  values = vapply(
    n,
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
