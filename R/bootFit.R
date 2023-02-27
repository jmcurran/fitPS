bootFit = function(x, B = 1000){
  yvals = rep(x$data$n, x$data$rn)
  n = length(yvals)

  to.psData = function(y){
    tbl = table(y)

    r = x
    r$n =  as.numeric(names(y))
    r$rn = tbl

    return(r)
  }

  boot.y = matrix(sample(yvals, n * B, replace = TRUE), nrow = B)
  boot.y = apply(boot.y, 1, to.psData)
  results = sapply(boot.y, function(y){
    fitDist(y)$shape
  })

  return(results)
}
