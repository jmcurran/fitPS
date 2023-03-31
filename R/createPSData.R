#' Manually create psData
#'
#' A mechanism for manually creating P or S data sets for use with fitDist.
#'
#' @param n a vector of labels for the number of groups or size of groups of
#'   glass, or a vector of observations.
#' @param rn if not \code{NULL} then this must be a vector of counts
#'   corresponding to each element of \code{n}. All entries must be greater than
#'   zero.
#' @param type either \code{"P"} or \code{"S"}.
#' @param notes a character string or a \code{\link[utils]{bibentry}}.
#'
#' @return an object of class \code{psData}. See \code{\link{readData}} for more
#'   details
#'
#' @examples
#' p = createPSData(0:2, c(98, 1, 1), type = "P")
#' p
#' @export
createPSData = function(n, rn = NULL, type = c("P", "S"), notes = NULL){

  type = match.arg(type)

  if(is.null(rn)){
    rn = table(n)
    n = as.numeric(names(rn))
    rn = as.numeric(rn)
  }else{
    if(length(n) != length(rn)){
      stop("n and rn must have equal length.")
    }
  }

  x = list(
    data = data.frame(n = n, rn = rn),
    type = type,
    notes = notes
  )

  class(x) = "psData"
  return(x)
}
