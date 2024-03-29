#' S3 method for objects of class \code{psData}
#'
#' Tests to see if two objects of class \code{psData} are equal. That is
#' their \code{type} is the same, and the data contained in \code{data} is the
#' same. See \code{\link{readData}} for a description of the \code{psData} class.
#'
#' @param lhs an object of class \code{psData}.
#' @param rhs an object of class \code{psData}.
#'
#' @details NOTE: the \code{notes} member variable is ignored in this function
#' as it is unlikely that a user would want to see if the notes are the same.
#'
#' @return TRUE if the two objects are equal
#' @export
#'
#' @examples
#' p = readData(system.file("extdata", "p.xlsx", package = "fitPS"))
#' p1 = makePSData(n = 0:2, count = c(98, 1, 1), type = "P")
#' p2 = makePSData(n = 0:2, count = c(97, 2, 1), type = "P")
#' p == p1 ## TRUE
#' p == p2 ## FALSE
#' p1 == p2 ## FALSE
`==.psData` = function(lhs, rhs){
  (lhs$type == rhs$type) && all(lhs$data == rhs$data)
}
