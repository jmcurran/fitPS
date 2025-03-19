fitdistLaplace = function(x, nterms = 10, ...){
  ## This is just proof of concept at the moment so
  ## the prior and its parameters are fixed to U[-2,2]
  zlpp = function(shape){
    -zetaloglikelihood(x, shape = shape) - log(zetaunifprior(shape))
  }

  start = 1.5

  fit = optim(par = start,
             fn = zlpp,
             method = "L-BFGS-B",
             lower = 0,
             hessian = TRUE)

  shape = fit$par ## NOTE VGAM's dzeta is parameterised in terms of
  ## s = alpha - 1. This has consequences in the formula below
  N = sum(x$data$rn)

  vshape = function(shape){
    z = VGAM::zeta(shape + 1)
    zprime = VGAM::zeta(shape + 1, 1)
    zprimeprime = VGAM::zeta(shape + 1, 2)
    numer = ((shape + 1) * z)^2
    denom = (shape + 1)^2 * N *(z * zprimeprime - zprime^2) + 1

    return(numer / denom)
  }

  var.shape = vshape(shape)

  fl = function(s){
    exp(zlpp(shape = shape)) *
      exp(-0.5 * (s - shape)^2 / var.shape)
  }

  return(fl)
}
