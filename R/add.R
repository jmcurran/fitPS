#' Add data to a psData object
#'
#' Add one or more new observations to an existing clothing survey object.
#'
#' @param x an object of class \code{psData}---see \code{\link{readData}} for details.
#' @param newData either a \code{vector}, \code{matrix} or \code{data.frame}
#' containing the new data. If a \code{vector} or \code{magtrix} is supplied then it must
#' be either of length or have two columns. If a \code{data.frame} is supplied then the columns
#' must be labelled \code{"n"} and \code{"rn"}. The new data MUST NOT contain values that already
#' exist in \code{x$n}
#'
#' @return an object of class \code{pSData}
#' @export
#'
#' @examples
#' add(Ssurveys$lau, c(11, 1))
add = function(x, newData){
  if(!is(x, "psData")){
    stop("x must be an object of class psData.")
  }

  if(!is(newData, "data.frame")){
    if(!is(newData, "numeric") && !is(newData, "matrix")){
      stop("newData must be a vector, a matrix, or a data.frame.")
    }
    if(is(newData, "numeric")){
      newData = data.frame(n = newData[1], rn = newData[2])
    }else{
      newData = as.data.frame(newData)
      names(newData) = c("n", "rn")
    }
  }

  if(ncol(newData) != 2){
    stop("newData must have exactly two columns.")
  }

  if(nrow(newData) == 0){
    stop("newData cannot be empty.")
  }

  ## make sure that the values of n are not overwriting existing data

  if(any(newData$n %in% x$data$n)){
    stop("Elements of newData already exist in the survey object.\nThis function does not support overwriting at this time.")
  }

  x$data = rbind(x$data, newData)
  return(x)
}
