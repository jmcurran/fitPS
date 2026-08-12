#' Validate that a zeta shape parameter is finite and greater than one.
#'
#' @param shape Zeta shape parameter on the fitPS scale.
#' @param argument Argument name used in validation error messages.
#' @return The validated shape value, invisibly or unchanged as used by callers.
#' @keywords internal
#' @noRd
validateZetaShape = function(shape, argument = "shape"){
  if(any(!is.finite(shape))){
    stop(argument, " must be finite.")
  }

  if(any(shape <= 1)){
    stop(argument, " must be greater than 1.")
  }

  invisible(shape)
}

#' Convert the fitPS zeta shape convention to the VGAM parameter convention.
#'
#' @param shape Zeta shape parameter on the fitPS scale.
#' @return The equivalent VGAM zeta parameter.
#' @keywords internal
#' @noRd
standardToVgamShape = function(shape){
  validateZetaShape(shape)
  shape - 1
}

#' Evaluate the zeta mass function using the fitPS shape convention.
#'
#' @param x An input object or numeric vector required by the helper.
#' @param shape Zeta shape parameter on the fitPS scale.
#' @param log Logical; return log probabilities when `TRUE`.
#' @return Numeric zeta probabilities or log probabilities.
#' @keywords internal
#' @noRd
dzetaStandard = function(x, shape, log = FALSE){
  vgamShape = standardToVgamShape(shape)
  VGAM::dzeta(x, shape = vgamShape, log = log)
}

#' Generate zeta random variates using the fitPS shape convention.
#'
#' @param n Number of random variates or observations.
#' @param shape Zeta shape parameter on the fitPS scale.
#' @return Integer-valued zeta random variates.
#' @keywords internal
#' @noRd
rzetaStandard = function(n, shape){
  vgamShape = standardToVgamShape(shape)
  VGAM::rzeta(n = n, shape = vgamShape)
}
