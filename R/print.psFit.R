#' S3 print method for an object of class \code{psFit}
#'
#' @param x an object of class \code{psFit}, usually from \code{\link{fitDist}}
#' @param ... other arguments passed to \code{print}
#'
#' @return No return value, called for side effects
#' @export
print.psFit = function(x, ...){
  cat(paste("The estimated shape parameter is ", round(x$shape + 1, 4), "\n"))
  cat(paste("The the standard error of shape parameter is ", round(sqrt(x$var.shape), 4), "\n"))
  cat(paste("The first ", length(x$fitted), "fitted values are:\n"))
  print(x$fitted)
}
