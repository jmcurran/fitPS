#' Generate zero inflated zeta random variates
#'
#' @param n the number of observations.
#' @param pi the mixing parameter for the zero-inflated zeta model---must be in
#'   (0, 1).
#' @param shape the shape parameter for the zero-inflated zeta. Must be greater
#'   than zero.
#'
#' @return a vector of random variates from a zero-inflated zeta model
#'
#' @examples
#' data(Psurveys)
#' roux = Psurveys$roux
#' fit.zi = fitZIDist(roux)
#' x = rZIzeta(n = sum(roux$data$rn), pi = fit.zi$pi, shape = fit.zi$shape)
#' table(x)
#' @export
rZIzeta = function(n, pi = 0.5, shape = 1){

  if(length(pi) > 1 || length(shape) > 1){
    stop("This function does not currently support vector valued inputs for Pi or shape.")
  }

  n = round(n)

  if(n <= 0){
    stop("n must be greater than zero.")
  }

  if(pi <= 0 || pi >= 1){
    stop("Pi must be in (0, 1)")
  }

  if(shape <= 0){
    stop("shape must be greater than zero.")
  }

  p = runif(n)
  x = p
  x[p < pi] = 1
  n = sum(p >= pi)
  x[p >= pi] = VGAM::rzeta(n, shape)
  return(x)
}
