#' Probability Functions
#'
#' Creates a probability function that allows the computation of any P or S
#' term.
#'
#' @param psFitobj an object of class \code{psFit}--see \code{\link{fitDist}}
#' and \code{\link{fitZIDist}}.
#'
#' @return a function that can be used to calculate any P or S term.
#' @export
#'
#' @examples
#' p = readData(system.file("extdata", "p.xlsx", package = "fitPS"))
#' fit = fitDist(p)
#' P = probfun(fit)
#' P(0:5)
probfun = function(psFitobj){
  pf = function(x){
    if(!psFitobj$zeroInflated){
      if(psFitobj$psData$type == "P"){
        p = dzeta(x + 1, shape = psFitobj$shape)
        names(p) = paste0("P", x)
        return(p)
      }else{
        p = dzeta(x, shape = psFitobj$shape)
        names(p) = paste0("S", x)
        return(p)
      }
    }else{
      if(psFitobj$psData$type == "P"){
        p = (1 - psFitobj$pi) * dzeta(x + 1, shape = psFitobj$shape)
        p[x == 0] = p[x == 0] + psFitobj$pi
        names(p) = paste0("P", x)
        return(p)
      }else{
        p = (1 - psFitobj$pi) * dzeta(x, shape = psFitobj$shape)
        p[x == 1] = p[x == 1] + psFitobj$pi
        names(p) = paste0("S", x)
        return(p)
      }
    }
  }
  return(pf)
}
