createPSData = function(n, rn, type = c("P", "S"), notes = NULL){

  type = match.arg(type)

  x = list(
    data = list(n = n, rn = rn),
    type = type,
    notes = notes
  )

  class(x) = "psData"
  return(x)
}
