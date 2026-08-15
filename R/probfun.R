#' Plug-in probability functions
#'
#' Creates a probability function that evaluates P or S terms through the
#' model contract at the fitted parameter values stored in a `psFit` object.
#'
#' @param psFitobj An object of class `psFit`.
#' @return A function that can calculate plug-in P or S terms.
#' @export
#'
#' @examples
#' p = readData(system.file("extdata", "p.xlsx", package = "fitPS"))
#' fit = fitDist(p)
#' P = probfun(fit)
#' P(0:5)
probfun = function(psFitobj) {
  if (!is(psFitobj, "psFit")) {
    stop("psFitobj must be an object of class psFit")
  }

  model = modelFromFit(psFitobj)
  parameters = fitModelParameters(psFitobj, model)
  surveyType = psFitobj$psData$type

  function(x) {
    probabilities = modelProbabilities(
      model = model,
      parameters = as.list(parameters),
      n = x,
      type = surveyType
    )
    result = as.numeric(probabilities[1L, ])
    names(result) = colnames(probabilities)
    result
  }
}
