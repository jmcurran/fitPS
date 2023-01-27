#' S3 method for objects of class \code{psData}
#'
#' Tests to see if two objects of class \code{psData} are equal. That is
#' their \code{type} is the same, and the data contained in \code{data} is the
#' same. See \code{\link{readData}} for a description of the \code{psData} class.
#'
#' @param lhs an object of class \code{psData}.
#' @param rhs an object of class \code{psData}.
#'
#' @return TRUE if the two objects are equal
#' @export
#'
#' @examples
#' p = readData(system.file("extdata", "p.xlsx", package = "fitPS"))
#' p1 = makePSData(n = 0:2, count = c(98, 1, 1), type = "P")
#' p == p1
`==.psData` = function(lhs, rhs){
  (lhs$type == rhs$type) && (lhs$data == rhs$data)
}
