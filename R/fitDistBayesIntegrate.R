fitDistBayesIntegrate = function(x, prior = makePrior(), nterms, ...){

  nvals = 1:nterms
  if(!is(x, "psData")){
    stop("x must be an object of class psData")
  }

  if(!is(x, "psData")){
    stop("x must be an object of class psData")
  }

  ## check which arguments have been provided and provide defaults
  ## if they haven't

  dotargs = list(...)
  is.arg = function(arg){
    return(tolower(arg) %in% tolower(names(dotargs)))
  }

  a = prior$range[1]
  b = prior$range[2]

  if(b <= a){
    stop("b must be greater than a!")
  }

  if(length(x$data$n) < 2){
    if(x$type == "S"){
      stop("There has to be at least one value higher than 1")
    }else{
      stop("There has to be at least one value higher than 0")
    }
  }

  obsData = if(x$type == 'P'){ ## the main difference is that the values need 1 added
    x$data$n + 1
  }else{
    x$data$n
  }

  logLik = function(shape){
    sum(x$data$rn * VGAM::dzeta(obsData, shape = shape, log = TRUE))
  }

  negLogLikTimesPrior = function(shape){
    - (logLik(shape) + prior$logd(shape))
  }

  # find the maximum of log likelihood times prior to determine a scaling factor
  # we scale with respect to the posterior mode
  opt <- optim(par = 1, fn = negLogLikTimesPrior, method = "L-BFGS-B")
  k0 <- opt$value

  s_posterior_unnormalised <- Vectorize(function(s)
    exp(-negLogLikTimesPrior(s) + k0))

  # integrate to determine the constant of proportionality K
  s_posterior_unnormalised_int <- integrate(s_posterior_unnormalised,
                                            lower = a, upper = b)
  logK <- log(s_posterior_unnormalised_int$value)

  # obtain true posterior such that the area under the curve is 1
  s_posterior <- Vectorize(function(s)
    exp(-negLogLikTimesPrior(s) + k0 - logK))

  # posterior mean and var by numerical integration
  s_posterior_mean_int <- integrate(f = function(s) s * s_posterior(s), lower = a, upper = b)
  s_posterior_mean2_int <- integrate(f = function(s) s^2 * s_posterior(s), lower = a, upper = b)

  s_posterior_variance <- s_posterior_mean2_int$value - s_posterior_mean_int$value^2

  fit = list()

  shape = s_posterior_mean_int$value
  var.shape = s_posterior_variance

  fit$par = shape

  fitted = VGAM::dzeta(nvals, shape = shape)
  names(fitted) = if(x$type == 'P'){
    paste0("P", nvals - 1)
  }else{
    paste0("S", nvals)
  }

  result = list(
    psData = x,
    fit = fit,
    shape = shape,
    var.shape = var.shape,
    fitted = fitted,
    model = "zeta",
    method = "integrate"
  )

  class(result) = "psFit"

  return(result)
}
