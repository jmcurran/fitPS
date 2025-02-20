#' S3 print method for an object of class \code{psFit}
#'
#' @param x an object of class \code{psFit}, usually from \code{\link{fitDist}}
#' or #' \code{\link{fitZIDist}}.
#' @param ... other arguments passed to \code{print}.
#'
#' @return No return value, called for side effects.
#' @export
print.psFit = function(x, ...){

  isBayes = x$method == "bayes"

  if(x$model == "ziz"){
    cat(paste("The estimated mixing parameter, pi, is", signif(x$pi, 4), "\n"))
  }

  if(x$model %in% c("zeta", "ziz")){
    if(isBayes){
      cat(paste("The estimated posterior mean of shape parameter is", round(x$shape + 1, 4), "\n"))
    }else{
      cat(paste("The estimated shape parameter is", round(x$shape + 1, 4), "\n"))
    }
  }

  if(x$model == "log"){
    cat(paste("The estimated model parameter pi is", signif(x$pi, 4), "\n"))
  }

  if(x$model == "zeta"){
    if(isBayes){
      cat(paste("The estimated posterior standard error of shape parameter is", round(sqrt(x$var.shape), 4), "\n"))
    }else{
      cat(paste("The standard error of shape parameter is", round(sqrt(x$var.shape), 4), "\n"))
    }
  }

  if(x$model %in% c("zeta", "ziz")){
    cat("------\n")
    writeLines(strwrap("NOTE: The shape parameter is reported so that it is consistent with Coulson et al. However, the value returned is actually s' = shape - 1 to be consistent with the VGAM parameterisation, which is used for computation. This has flow on effects, for example in confInt. This will be changed at some point.\n"))
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
}
