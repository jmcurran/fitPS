#' Generate random variates from a zeta distribution
#'
#' @param n Same as \code{\link[stats]{Poisson}}.
#' @param shape The standard zeta shape parameter, greater than 1.
#'
#' See \code{\link[VGAM]{rzeta}}.
#'
#' @export
rzeta = function(n, shape){
  rzetaStandard(n = n, shape = shape)
}
