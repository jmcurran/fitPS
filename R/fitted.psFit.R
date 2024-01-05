#' S3 fitted method for an object of class \code{psFit}
#'
#' @param object an object of class \code{psFit}, usually from \code{\link{fitDist}}
#' or \code{\link{fitZIDist}}.
#' @param n This parameter is \code{NULL} by default. If it is not \code{NULL}
#'  then it must be either the number of fitted terms to be return, or, a vector
#' containing the desired fitted values.
#' @param ... other arguments passed to \code{fitted}---not used.
#'
#' @return a named vector of fitted probabilities
#' @importFrom stats fitted
#' @export
fitted.psFit = function(object, n = NULL, ...){
  if(is.null(n)){
    return(object$fitted)
  }else{
    if(any(n < 0)){
      stop("The value(s) of n must be >= 0")
    }

    if(any(abs(n - floor(n) > .Machine$double.eps))){
      n = floor(n)
      warning("n must be integer valued, and so has been rounded down")
    }

    type = object$psData$type

    if(type == "S" && any(n == 0)){
      stop("n must be >= 1 for S terms")
    }

    pf = probfun(object)

    if(length(n) == 1){
      start = ifelse(type == "P", 0, 1)
      end = n - ifelse(type == "P", 1, 0)
      rval = pf(start:end)
    }else{
      rval = pf(sort(n))
    }

    return(rval)
  }
}
