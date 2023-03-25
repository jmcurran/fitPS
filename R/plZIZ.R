#' @importFrom ellipse ellipse
plZI  = function(x, level){
  if(!x$zeroInflated){
    stop("This function is for the ZIZ model only")
  }

  v = fit$var.cov


}
