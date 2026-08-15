#' Converts an object of class \code{psData} to a \code{data.frame}
#'
#' Converts an object of class \code{psData}---see \code{\link{readData}}---to a
#' \code{data.frame} that can be used with in functions in other packages such as
#' \code{\link[VGAM]{vglm}} to fit more complicated models.
#'
#' @param x an object of class \code{psData}---see \code{\link{readData}}for more
#' details.
#' @param ... any other arguments passed to \code{data.frame}.
#'
#' @details If \code{x} is a \code{psData} object of type \code{"P"}, i.e. it
#'   relates to numbers of groups of glass, then a \code{data.frame} with a single variable
#'   \code{count} will be return where \code{count = rep(x$data$n + 1,
#'   x$data$rn)}. The counts have one added to them because the zeta
#'   distribution requires that the counts are greater than or equal to one.  If
#'   \code{x} is a \code{psData} object of type \code{"P"}, i.e. it relates to
#'   group sizes, then a \code{data.frame} with a single variable \code{count}
#'   will be return where \code{count = rep(x$data$n, x$data$rn)}.
#'
#' @return a \code{data.frame} with a single variable \code{count}. The number
#' of rows in the \code{data.frame} is equal to \code{sum(x$data$rn)}.
#'
#' @export
#'
#' @examples
#' p = readData(system.file("extdata", "p.xlsx", package = "fitPS"))
#' p.df = as.data.frame(p)
#' table(p.df$count)
#' p$data
as.data.frame.psData = function(x, ...){
  if(x$type == "P"){
    data.frame(count = rep(x$data$n + 1, x$data$rn), ...)
  }else{
    data.frame(count = rep(x$data$n, x$data$rn), ...)
  }
}
#' Variance generic
#'
#' @param x an object for which we want to compute the sample variance.
#' @param \dots Any additional arguments to be passed to \code{var}.
#' @export
var = function(x, ...){
  UseMethod("var")
}
#' @describeIn var Default variance method for non-psData objects.
#' @export
var.default = function(x, ...){
  stats::var(x, ...)
}
#' An S3 method for computing the mean of clothing survey for the number of
#' groups or size of groups
#'
#' @param x an object of class \code{psData}---\code{\link{readData}} for more
#'   details.
#' @param ... other arguments which are passed to \code{\link[base]{sum}}
#'
#' @return the mean of the data. If there are \eqn{r_i}{r[i]} observations of
#'   the value \eqn{n_i}{n[i]} then the mean is given by
#'   \deqn{\sum_i\frac{r_i\times n_i}{\sum_i{r_i}}}{sum(r[i]*n[i])/sum(r[i])}.
#' @export
#'
#' @examples
#' data(Psurveys)
#' mean(Psurveys$roux)
mean.psData = function(x, ...){
  return(sum(x$data$n * x$data$rn, ...) / sum(x$data$rn, ...))
}
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
#' S3 print method for an object of class \code{psData}
#'
#' @param x an object of class \code{psData}, usually from \code{\link{readData}}
#' or \code{\link{makePSData}}
#' @param ... other arguments passed to \code{print}
#'
#' @importFrom knitr kable
#' @importFrom utils capture.output

#' @return No return value, called for side effects
#' @export
print.psData = function(x, ...){
  kbl = kable(x$data, format = "simple",
              label = NA,
              caption = ifelse(x$type == "P",
              "Number of Groups",
              "Group Size")
             )

  kbl[1] = gsub("^Table[:] +(.*$)", "\\1", kbl[1])

  z = capture.output(print(kbl, ...)) ## There might be a smarter way to get rid of kables two LF characters at the start of print, but I don't know what it is.
  z = z[-c(1:2)]
  cat(z, sep = "\n")

  if(!is.null(x$notes)){
    print(x$notes, ...)
  }
}
#' S3 method for objects of class \code{psData}
#'
#' Tests to see if two objects of class \code{psData} are equal. That is
#' their \code{type} is the same, and the data contained in \code{data} is the
#' same. See \code{\link{readData}} for a description of the \code{psData} class.
#'
#' @param lhs an object of class \code{psData}.
#' @param rhs an object of class \code{psData}.
#'
#' @details NOTE: the \code{notes} member variable is ignored in this function
#' as it is unlikely that a user would want to see if the notes are the same.
#'
#' @return TRUE if the two objects are equal
#' @export
#'
#' @examples
#' p = readData(system.file("extdata", "p.xlsx", package = "fitPS"))
#' p1 = makePSData(n = 0:2, count = c(98, 1, 1), type = "P")
#' p2 = makePSData(n = 0:2, count = c(97, 2, 1), type = "P")
#' p == p1 ## TRUE
#' p == p2 ## FALSE
#' p1 == p2 ## FALSE
`==.psData` = function(lhs, rhs){
  (lhs$type == rhs$type) && all(lhs$data == rhs$data)
}
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
