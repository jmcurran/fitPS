#' Create a survey data set manually
#'
#' Create a survey data set from the command line rather than reading data in
#' from a file. This function is likely to be only useful where there are a very
#' small number of group sizes, or sizes of groups of glass.
#'
#' @param n Either the number of groups of glass or the size of different groups
#' of glass.
#' @param count Either the number of people in the survey sample who had \eqn{n}{n}
#' groups of glass on their clothing, or the number of people who had a group of glass of
#' size \eqn{n}{n}.
#' @param type either \code{"P"} or \code{"S"}
#'
#' @return an object of type \code{psData}---see \code{\link{readData}} for more
#' details.
#'
#' @aliases makeData
#' @export
#'
#' @examples
#' ## recreate the data read in the readData example
#' p1 = makePSData(n = c(0, 1, 2), count = c(98, 1, 1), type = "P")
#' s1 = makePSData(n = 1:3, count = c(1, 1, 1), type = "S")
#' p1
#' s1
#'
#' @seealso readData
makePSData = makeData = function(n, count, type = c("P", "S")){
  type = match.arg(type)

  dataf = data.frame(n = n, rn = count)

  if(any(dataf$rn <= 0)){
    stop("You must have at least one non-zero positive count")
  }

  if(type == "S" & any(dataf$n <= 0)){
    stop("The labels for group sizes must be 1 or greater")
  }else if(type == "P" & any(dataf$n < 0)){
    stop("The labels for group sizes must be 1 or greater")
  }

  if(any(dataf$n -floor(dataf$n) > 0) || any(dataf$rn -floor(dataf$rn) > 0)){
    dataf$n = round(dataf$n)
    dataf$rn = round(dataf$rn)
    warning("Non-integer input detected. Both input vectors have been rounded.")
  }

  result = list(type = type,
                data = dataf)
  class(result) = "psData"

  return(result)

}
