#' Variance generic
#'
#' @param x an object for which we want to compute the sample variance.
#' @param \dots Any additional arguments to be passed to \code{var}.
#' @export
var = function(x, ...){
  UseMethod("var")
}
#' @export
var.default = function(x, ...){
  stats::var(x, ...)
}
#' An S3 method for computing the mean of clothing survey for the number of
#' groups or size of groups
#'
#' @param x an object of class \code{psData}---\code{\link{readData}} for more
#'   details.
#' @param ... other arguments which are passed to \code{\link[base]{sum}}
#'
#' @return the mean of the data. If there are \eqn{r_i}{r[i]} observations of
#'   the value \eqn{n_i}{n[i]} then the mean is given by
#'   \deqn{\sum_i\frac{r_i\times n_i}{\sum_i{r_i}}}{sum(r[i]*n[i])/sum(r[i])}.
#' @export
#'
#' @examples
#' data(Psurveys)
#' mean(Psurveys$roux)
mean.psData = function(x, ...){
  return(sum(x$data$n * x$data$rn, ...) / sum(x$data$rn, ...))
}
