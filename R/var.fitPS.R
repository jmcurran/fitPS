#' An S3 method for computing the variance of clothing survey for the number of
#' groups or size of groups
#'
#' @param x an object of class \code{psData}---\code{\link{readData}} for more
#'   details.
#' @param ... other arguments which are passed to \code{\link[base]{sum}}
#'
#' @return the mean of the data. If there are \eqn{r_i}{r[i]} observations of
#'   the value \eqn{n_i}{n[i]} then the variance is computed by
#'   \eqn{\mathrm{E}[X^2]-\mathrm{E}[X]^2}{E[X^2]-E[X]^2}, where
#'   \eqn{\mathrm{E}[X]}{E[X]} is computed using \deqn{\sum_i\frac{r_i\times
#'   n_i}{\sum_i{r_i}}}{sum(r[i]*n[i])/sum(r[i])} , and
#'   \eqn{\mathrm{E}[X^2]}{E[X^2]} is computed by \deqn{\sum_i\frac{r_i\times
#'   n_i^2}{\sum_i{r_i}}}{sum(r[i]*n[i]^2)/sum(r[i])}. We realise that the
#'   computational formula,
#'   \eqn{\mathrm{E}[X^2]-\mathrm{E}[X]^2}{E[X^2]-E[X]^2}, is usually not
#'   regarded as computationally stable, but the magnitude of the numbers
#'   involved is such that, that this is not likely to cause an issue.
#'
#' @examples
#' data(Psurveys)
#' var(Psurveys$roux)
#' @export
var.psData = function(x, ...){
  Ex.sq = mean(x, ...)^2
  Ex2 = (sum(x$data$n^2 * x$data$rn, ...) / sum(x$data$rn, ...))
  return(Ex2 - Ex.sq)
}
