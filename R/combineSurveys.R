combineSurveys = function(...){

  Surveys = list(...)
  n = lapply(Surveys, function(x){
    rep(x$data$n, x$data$rn)
    })

  return(makePSData(unlist(n)))
}
