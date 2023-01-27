#' Converts an object of class \code{psData} to a \code{data.frame}
#'
#' Converts an object of class \code{psData}---see \code{\link{readData}}---to a
#' \code{data.frame} that can be used with \code{VGAM}} to fit more complicated
#' models.
#'
#' @param x an object of class \code{psData}---see \code{\link{readData}}for more
#' details.
#' @param ... any other arguments passed to \code{data.frame}.
#'
#' @details If \code{x} is a \code{psData} object of type \code{"P"}, i.e. it
#'   relates to numbers of groups of glass, then a \code{data.frame} with a single variable
#'   \code{count} will be return where \code{count = rep(x$data$n + 1,
#'   x$data$rn)}. The counts have one added to them because the zeta
#'   distribution requires that the counts are greater than or equal to one.  If
#'   \code{x} is a \code{psData} object of type \code{"P"}, i.e. it relates to
#'   group sizes, then a \code{data.frame} with a single variable \code{count}
#'   will be return where \code{count = rep(x$data$n, x$data$rn)}.
#'
#' @return a \code{data.frame} with a single variable \code{count}. The number
#' of rows in the \code{data.frame} is equal to \code{sum(x$data$rn)}.
#'
#' @export
#'
#' @examples
#' p = readData(system.file("extdata", "p.xlsx", package = "fitPS"))
#' p.df = as.data.frame(p)
#' table(p.df$count)
#' p$data
as.data.frame.psData = function(x, ...){
  if(x$type == "P"){
    data.frame(count = rep(x$data$n + 1, x$data$rn), ...)
  }else{
    data.frame(count = rep(x$data$n, x$data$rn), ...)
  }
}
