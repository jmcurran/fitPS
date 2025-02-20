fitDistBayes = function(x, nterms, ...){

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
  a = ifelse(is.arg("a"), dotargs$a, -2)
  b = ifelse(is.arg("b"), dotargs$b, 2)
  nIter = ifelse(is.arg("nIter"), dotargs$nIter, 1e4)
  nBurnIn = ifelse(is.arg("nBurnIn"), dotargs$nBurnIn, 1e3)
  silent = ifelse(is.arg("silent"), dotargs$silent, TRUE)

  if(nIter < 1000){
    warning("The number of samples from the MCMC chain really should be 1000 or higher.")
  }

  if(nIter <= 0 || nBurnIn <=0){
    stop("nIter and nBurnIn must be greater than zero.")
  }

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

  log.shape = runif(nTotal, a, b)
  shape0 = max(shape0 - 1, .Machine$double.eps)
  draws = exp(-log.shape)
  shape1 = draws[1]

  chain = rep(shape0, nIter)
  log.u = log(runif(nTotal))

  ll0 = logLik(shape0) + log(W * shape0)
  i = 1

  if(!silent){
    pb = txtProgressBar(1, nTotal, 1, style = 3, label = 'Burning in')
  }

  while(i <= nTotal){
    if(i <= nTotal){
      shape1 = draws[i]
      ll1 = logLik(shape1) + log(W * shape1)
    }

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
