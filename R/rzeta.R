#' Generate random variates from a zeta distribution
#'
#' @param n Same as \code{\link[stats]{Poisson}}.
#' @param shape The positive shape parameter \deqn{s = \alpha - 1}{s = alpha - 1}.
#'
#' See \code{\link[VGAM]{rzeta}}.
#'
#' @export
rzeta = function(n, shape){
  VGAM::rzeta(n = n, shape = shape)
}
