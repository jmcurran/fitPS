#' Read count data from file
#'
#' Reads observed counts of either the number of groups or the size of the
#' groups. The file must have only two columns. One of the columns must be
#' labelled P or S and the other count. It does not matter if the column names
#' are in upper case or not. The P column can have labels 0, 1, 2, \ldots
#' representing the observation of 0, 1, 2, or more groups. The corresponding
#' count column should contain a positive (non-zero) count for each number of
#' groups. Similarly, if the file contains S counts, then the S column can
#' contain labels 1, 2, \ldots representing the observation of 1, 2, \ldots
#' fragments in a group. Note that zeros are neither allowed, or useful, in the
#' file as they both simply result in log-likelihood terms of zero, and
#' therefore make no difference.
#'
#' @param fileName the name of the file to be read. Must be either a modern
#'   (xlsx) Excel file or a csv file.
#' @param notes any additional information about the data, such as the source or
#'   a reference.
#' @param ... any additional parameters which will be passed to either
#'   \code{read_excel} or \code{read.csv} depending on the extension of your
#'   input file.
#'
#' @return an object of class \code{psData} which is a list containing member
#'   variables:
#' \itemize{
#'   \item{\code{type}}{ -- either \code{"P"} or \code{"S"}}
#'   \item{\code{data}}{ -- a \code{data.frame} which contains columns
#' \code{n} and \code{rn}, representing the number of groups/fragments, and the
#' number of times that was seen, respectively.}
#'   \item{\code{notes}}{ --- either a \code{\link[utils]{bibentry}} or a character
#'   string which allows extra information about the data to be stored, such as
#'   the source, or reference.}
#' }
#'
#' @import dplyr
#' @importFrom utils read.csv
#' @importFrom readxl read_excel
#'
#' @export
#'
#' @examples
#' p = readData(system.file("extdata", "p.xlsx", package = "fitPS"))
#' p
#' s = readData(system.file("extdata", "s.xlsx", package = "fitPS"))
#' s
readData = function(fileName, notes = NULL, ...){
  extn = tools::file_ext(fileName)

  if(extn == "xlsx"){
    dataf = readxl::read_excel(fileName, ...)
  }else if(extn == "csv"){
    dataf = read.csv(fileName, ...)
  }else{
    stop("I don't know how to read this data")
  }

  if(ncol(dataf) != 2){
    stop("I'm expecting only two columns")
  }

  vars = tolower(names(dataf))

  if(!("count" %in% vars)){
    stop("You need a column named count in your data")
  }

  vars = vars[vars != "count"]

  if(!(tolower(vars) %in% c("p", "s"))){
    stop("One of your columns should be labelled P or S")
  }

  if(nrow(dataf) < 1){
    stop("You need at least one observation to fit the distribution")
  }

  result = data.frame(rn = as.integer(dataf$count))

  result = result |>
    mutate(n = dataf |>
            select(-matches("^count$")) |>
             pull() |>
             as.integer()
             )
  # result = result |>
  #   dplyr::select(c(n, rn)) ## rearrange the order of the variables
  ## The code above works, but somehow bothers R CHECK. Going for the blunt
  ## approach
  result = result[, c(2, 1)]

  type = if(vars == "p"){
    "P"
  }else{
    "S"
  }

  if(length(unique(result$n)) < length(result$n)){
    stop("The categories (n values) must be unique")
  }

  ## Final test - check all the counts are positive integers
  if(any(result$n < 0) || any(result$rn < 0)){
    stop("The categories and counts must be >= 0")
  }

  if(any(result$rn <= 0)){
    stop("You must have at least one non-zero positive count")
  }

  result = list(type = type,
                data = result,
                notes = notes)
  class(result) = "psData"

  return(result)
}
