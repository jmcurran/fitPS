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
    sum(x$data$rn * dzetaStandard(obsData, shape = shape, log = TRUE))
  }

  negLogLikTimesPrior = function(shape){
    - (logLik(shape) + prior$logd(shape))
  }

  # find the maximum of log likelihood times prior to determine a scaling factor
  # we scale with respect to the posterior mode
  opt = optim(
    par = mean(c(a, b)),
    fn = negLogLikTimesPrior,
    method = "L-BFGS-B",
    lower = a,
    upper = b
  )
  k0 = opt$value

  sPosteriorUnnormalised = Vectorize(function(s) {
    exp(-negLogLikTimesPrior(s) + k0)
  })

  # integrate to determine the constant of proportionality K
  sPosteriorUnnormalisedInt = integrate(
    sPosteriorUnnormalised,
    lower = a,
    upper = b
  )
  logK = log(sPosteriorUnnormalisedInt$value)

  # obtain true posterior such that the area under the curve is 1
  sPosterior = Vectorize(function(s) {
    exp(-negLogLikTimesPrior(s) + k0 - logK)
  })

  # posterior mean and var by numerical integration
  sPosteriorMeanInt = integrate(f = function(s) {
    s * sPosterior(s)
  }, lower = a, upper = b)
  sPosteriorMean2Int = integrate(f = function(s) {
    s^2 * sPosterior(s)
  }, lower = a, upper = b)

  sPosteriorVariance = sPosteriorMean2Int$value - sPosteriorMeanInt$value^2

  fit = list()

  shape = sPosteriorMeanInt$value
  var.shape = sPosteriorVariance

  fit$par = shape

  fitted = dzetaStandard(nvals, shape = shape)
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
    pdf = sPosterior,
    model = "zeta",
    method = "integrate"
  )

  class(result) = "psFit"

  return(result)
}
