bootFit = function(x, B = 1000, model = c("zeta", "zi.zeta")){
  yvals = rep(x$data$n, x$data$rn)
  n = length(yvals)

  to.psData = function(y){
    tbl = table(y)
    counts = as.vector(tbl)

    r = x
    r$data = data.frame(n =  as.numeric(names(tbl)),
                        rn =  counts)

    return(r)
  }

  boot.y = matrix(sample(yvals, n * B, replace = TRUE), nrow = B)
  boot.y = apply(boot.y, 1, to.psData)

  model = match.arg(model)

  if(model == "zeta"){
    results = sapply(boot.y, function(y){
      fitDist(y)$shape
    })
  }else{
    results = lapply(boot.y, function(y){
      fit = fitZIDist(y)
      return(c(fit$pi, fit$shape))
    })
    results = do.call("rbind", results)
  }

  return(results)
}
