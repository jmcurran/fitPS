validateZetaShape = function(shape, argument = "shape"){
  if(any(!is.finite(shape))){
    stop(argument, " must be finite.")
  }

  if(any(shape <= 1)){
    stop(argument, " must be greater than 1.")
  }

  invisible(shape)
}

standardToVgamShape = function(shape){
  validateZetaShape(shape)
  shape - 1
}

dzetaStandard = function(x, shape, log = FALSE){
  vgamShape = standardToVgamShape(shape)
  VGAM::dzeta(x, shape = vgamShape, log = log)
}

rzetaStandard = function(n, shape){
  vgamShape = standardToVgamShape(shape)
  VGAM::rzeta(n = n, shape = vgamShape)
}
