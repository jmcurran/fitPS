fitDistBayes = function(x, prior = makePrior(), nterms, ...){

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

  shape0 = ifelse(is.arg("shape0"), dotargs$shape0, 1)

  nIter = ifelse(is.arg("nIter"), dotargs$nIter, 1e4)
  nBurnIn = ifelse(is.arg("nBurnIn"), dotargs$nBurnIn, 1e3)
  silent = ifelse(is.arg("silent"), dotargs$silent, TRUE)

  if(nIter < 1000){
    warning("The number of samples from the MCMC chain really should be 1000 or higher.")
  }

  if(nIter <= 0 || nBurnIn <=0){
    stop("nIter and nBurnIn must be greater than zero.")
  }

  a = prior$range[1]
  b = prior$range[2]

  if(b <= a){
    stop("b must be greater than a!")
  }

  W = b - a
  if(W <= 0){
    stop("This should never happen.")
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
    #-sum(VGAM::dzeta(rep(obsData, x$data$rn), shape = shape, log = TRUE))
    sum(x$data$rn * VGAM::dzeta(obsData, shape = shape, log = TRUE))
  }

  nTotal = nIter + nBurnIn

  draws = runif(nTotal, a, b)

  chain = numeric(nIter)
  log.u = log(runif(nTotal))

  priorLogd <- prior$logd
  ll0 = logLik(shape0) + priorLogd(shape0)

  if (!is.finite(ll0)){
    stop("Log likelihood is not finite at starting value")
  }

  if(!silent){
    pb = txtProgressBar(1, nTotal, 1, style = 3, label = 'Burning in')
  }

  i = 1
  while(i <= nTotal){

    shape1 = draws[i]
    ll1 = logLik(shape1) + priorLogd(shape1)

    if(ll1 > ll0 || log.u[i] < (ll1 - ll0)){
      shape0 = shape1
      ll0 = ll1
    }

    if(i > nBurnIn){
      chain[i - nBurnIn] = shape0
    }
    i = i + 1
    if(!silent){
      if(i <= nBurnIn){
        setTxtProgressBar(pb, i)
      }else{
        setTxtProgressBar(pb, i, label = 'Sampling')
      }
    }
  }

  if(!silent){
    close(pb)
  }

  fit = list()
  fit$par = shape = mean(chain)
  var.shape = var(chain)

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
    chain = chain,
    model = "zeta",
    method = "bayes"
  )

  class(result) = "psFit"

  return(result)
}
