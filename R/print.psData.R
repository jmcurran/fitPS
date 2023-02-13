#' S3 print method for an object of class \code{psData}
#'
#' @param x an object of class \code{psData}, usually from \code{\link{readData}}
#' or \code{\link{makePSData}}
#' @param ... other arguments passed to \code{print}
#'
#' @importFrom knitr kable
#' @return No return value, called for side effects
#' @export

print.psData = function(x, ...){
  print(knitr::kable(x$data, format = "simple"), ...)
  if(x$type == "P"){
    cat("\nNumber of groups data\n\n")
  }else{
    cat("\nGroup size data\n\n")
  }
  if(!is.null(x$notes)){
    print(x$notes, ...)
  }
}
