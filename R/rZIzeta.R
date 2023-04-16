#' Generate zero inflated zeta random variates
#'
#' @param n the number of observations.
#' @param pi the mixing parameter for the zero-inflated zeta model---must be in
#'   (0, 1).
#' @param shape the shape parameter for the zero-inflated zeta. Must be greater
#'   than zero.
#' @param offset the zeta distribution returns random variates that are greater
#'   than, or equal to one. If the offset is greater than 0, then the
#'   distribution is anchored on (has minimum value of) \code{1 - offset}.
#'
#' @return a vector of random variates from a zero-inflated zeta model
#'
#' @details Technically this function returns values from the one-inflated zeta
#'   distribution. However, if \code{offset} is greater than zero (and typically
#'   we expect it to be 1), then the minimium random variate value is \code{1 -
#'   offset}. We chose the name "zero-inflated zeta" as more people are familiar
#'   with zero-inflated models.
#'
#'
#' @examples
#' data(Psurveys)
#' roux = Psurveys$roux
#' fit.zi = fitZIDist(roux)
#' x = rZIzeta(n = sum(roux$data$rn), pi = fit.zi$pi, shape = fit.zi$shape)
#' table(x)
#' @export
rZIzeta = function(n, pi = 0.5, shape = 1, offset = 0){

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
  return(x - offset)
}

#' @rdname rZIzeta
#' @export
rzizeta = rZIzeta

#' @rdname rZIzeta
#' @export
rzizeta = rZIzeta
rziz = rZIzeta
