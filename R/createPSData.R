#' Manually create psData
#'
#' A mechanism for manually creating P or S data sets for use with fitDist.
#'
#' @param n a vector of lables for the number of groups or size of groups of
#'   glass.
#' @param rn a vector of counts corresponding to each element of \code{n}. All
#'   entries must be greater than zero.
#' @param type either \code{"P"} or \code{"S"}.
#' @param notes a character string or a \code{\link[utils]{bibentry}}.
#'
#' @return an object of class \code{psData}. See \code{\link{readData}} for more
#'   details
#'
#' @examples
#' p = createData(0:2, c(98, 1, 1), type = "P")
#' p
#' @export
createPSData = function(n, rn, type = c("P", "S"), notes = NULL){

  type = match.arg(type)

  x = list(
    data = data.frame(n = n, rn = rn),
    type = type,
    notes = notes
  )

  class(x) = "psData"
  return(x)
}
