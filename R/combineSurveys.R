#' Combine compatible P/S survey objects into one data object.
#'
#' @param ... Additional arguments passed to the underlying fitting or helper routine.
#' @return A combined `psData` object.
#' @keywords internal
#' @noRd
combineSurveys = function(...){

  Surveys = list(...)
  n = lapply(Surveys, function(x){
    rep(x$data$n, x$data$rn)
    })

  return(makePSData(unlist(n)))
}
