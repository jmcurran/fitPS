#' S3 print method for an object of class \code{psFit}
#'
#' @param x an object of class \code{psFit}, usually from \code{\link{fitDist}}
#' @param ... other arguments passed to \code{print}
#'
#' @return No return value, called for side effects
#' @export
print.psFit = function(x, ...){
  cat(paste("The estimated shape parameter is ", round(x$shape + 1, 4), "\n"))
  cat(paste("The standard error of shape parameter is ", round(sqrt(x$var.shape), 4), "\n"))
  cat("------\n")
  cat("NOTE: The shape parameter is reported so that it is consistent with Coulson et al.\n")
  cat("However, the value returned is actually s' = shape - 1 to be consistent with the \n")
  cat("VGAM parameterisation, which is used for computation. This has flow on effects, for\n")
  cat("example in confInt. This will be changed at some point.\n")
  cat("------\n\n")

  args = list(...)
  if("nterms" %in% names(args)){
    nterms = as.integer(args$nterms[1])

    if(nterms < 1){
      stop("nterms must be >= 1")
    }else if(nterms > 10){
      nvals = 1:nterms
      fitted = VGAM::dzeta(nvals, shape = x$shape)
      names(fitted) = if(x$type == 'P'){
        paste0("P", nvals - 1)
      }else{
        paste0("S", nvals)
      }
      cat(paste("The first ", nterms, "fitted values are:\n"))
      print(fitted)
    }else{
      cat(paste("The first ", nterms, "fitted values are:\n"))
      print(x$fitted[1:nterms])
    }
  }else{
    cat(paste("The first ", length(x$fitted), "fitted values are:\n"))
    print(x$fitted)
  }
}
