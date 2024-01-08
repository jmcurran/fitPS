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
  kbl = knitr::kable(x$data, format = "simple",
                     label = NA,
                     caption = ifelse(x$type == "P",
                                "Number of Groups",
                                "Group Size")
                    )

  kbl[1] = gsub("^Table[:] +(.*$)", "\\1", kbl[1])

  z = capture.output(print(kbl, ...))
  z = z[-c(1:2)]
  cat(z, sep = "\n")

  if(!is.null(x$notes)){
    print(x$notes, ...)
  }
}
