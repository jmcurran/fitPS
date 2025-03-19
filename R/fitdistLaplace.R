fitdistLaplace = function(x, nterms = 10, s, ...){
  fit = fitdist(x, nterms = nterms, method = "mle", ...)
  shape.hat = fit$shape
  sigma.sq.hat = fit$var.shape

  fl = function(s){
    fitsPS:::zetalikelihood(x, shape = s) *
      exp(-0.5 * (s - shape.hat)^2 / sigma.sq.hat)
  }

  return(fl)
}
