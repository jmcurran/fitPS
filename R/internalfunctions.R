#' Evaluate the zeta log likelihood for internal Bayesian calculations.
#'
#' @param x An input object or numeric vector required by the helper.
#' @param shape Zeta shape parameter on the fitPS scale.
#' @return A numeric log-likelihood value.
#' @keywords internal
#' @noRd
zetaloglikelihood = function(x, shape){
  offset = ifelse(x$type == 'P', 1, 0)
  if(length(shape) > 1){
    ll = function(s){
      sum(x$data$rn * dzetaStandard(x$data$n + offset, shape = s, log = TRUE))
    }
    sapply(shape, ll)
  }else{
    sum(x$data$rn * dzetaStandard(x$data$n + offset, shape = shape, log = TRUE))
  }
}

zll = zetaloglikelihood

#' Evaluate the zeta likelihood on the natural probability scale.
#'
#' @param x An input object or numeric vector required by the helper.
#' @param shape Zeta shape parameter on the fitPS scale.
#' @return A numeric likelihood value.
#' @keywords internal
#' @noRd
zetalikelihood = function(x, shape){
  exp(zll(x, shape))
}

#' Evaluate the legacy bounded uniform prior density for the zeta shape.
#'
#' @param s Candidate zeta shape value.
#' @param a Lower bound of the legacy uniform prior.
#' @param b Upper bound of the legacy uniform prior.
#' @return A numeric prior density.
#' @keywords internal
#' @noRd
zetaunifprior = function(s, a = exp(-2), b = exp(2)){
  result = numeric(length(s))
  inRange = s >= a & s <= b
  result[inRange] = 1 / ((b - a) * s[inRange])
  result[!inRange] = 0
  return(result)
}

zup = zetaunifprior
