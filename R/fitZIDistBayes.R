
#' @importFrom stats cov rbeta dbeta
fitZIDistBayes = function(x, nterms = 10,
                      ...){
  nvals = 1:nterms
  if(!is(x, "psData")){
    stop("x must be an object of class psData")
  }

  ## check which arguments have been provided and provide defaults
  ## if they haven't

  dotargs = list(...)
  is.arg = function(arg){
    return(tolower(arg) %in% tolower(names(dotargs)))
  }

  theta0 = if(is.arg("theta0")){
             dotargs$theta0
           }else{
             c(0.5, 1)
           }
  a = ifelse(is.arg("a"), dotargs$a, -2)
  b = ifelse(is.arg("b"), dotargs$b, 2)
  shape1 = ifelse(is.arg("shape1"), dotargs$shape1, 1)
  shape2 = ifelse(is.arg("shape2"), dotargs$shape2, 1)
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

  if(theta0[1] <= 0 || theta0[1] >= 1){
    stop("The starting value for pi must be in (0, 1)")
  }

  if(theta0[2] <= 0){
    stop("The Zeta function is undefined for shape = 0. Choose a start value > 0.")
  }

  if(b <= a){
    stop("b must be greater than a!")
  }

  W = b - a
  if(W <= 0){
    stop("This should never happen.")
  }

  if(shape1 <= 0 || shape2 <= 0){
    stop("shape1, shape2 must be greater than zero.")
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

  y = rep(obsData, x$data$rn)

  d.one.inflated.zeta = function(x, shape, p, log = FALSE){
    rval = (1 - p) * VGAM::dzeta(x, shape = shape)
    rval[x == 1] = rval[x == 1] + p

    if(log){
      return(log(rval))
    }
    return(rval)
  }

  logLik = function(params){
    p = params[1]
    shape = params[2]

    rval = (1 - p) * VGAM::dzeta(obsData, shape = shape)
    rval[obsData == 1] = rval[obsData == 1] + p

    r = sum(x$data$rn * log(rval))
    if(is.infinite(r) || is.nan(r)){
      stop(sprintf("Infinite log-likelihod: pi = %6.4E shape = %6.4f\n", p, shape))
    }
    return(r)
  }

  nTotal = nIter + nBurnIn

  s0 = max(theta0[2] - 1, .Machine$double.eps)
  s.draws = runif(nTotal, exp(a), exp(b))
  s1 = s.draws[1]

  pi0 = min(max(theta0[1], .Machine$double.eps), 1 - .Machine$double.eps)
  pi.draws = if(shape1 == 1 && shape2 == 1){
              runif(nTotal)
             }else{
              rbeta(nTotal, shape1, shape2)
             }
  pi1 = pi.draws[1]

  whichParam = sample(1:2, nTotal, TRUE)

  chain = data.frame(pi = rep(pi0, nIter), shape = rep(s0, nIter))
  log.u = log(runif(nTotal))

  ll0 = logLik(c(pi0, s0)) - log(W * s0) + dbeta(pi0, shape1, shape2, log = TRUE)
  i = 1

  if(!silent){
    pb = txtProgressBar(1, nTotal, 1, style = 3, label = 'Burning in')
  }

  while(i <= nTotal){
    if(i <= nTotal){
      if(whichParam[i] == 1){
        s1 = s.draws[i]
      }else{
        pi1 = pi.draws[i]
      }
      ll1 = logLik(c(pi1, s1)) - log(W * s1) - dbeta(pi1, shape1, shape2, log = TRUE)
    }

    if(ll1 > ll0 || log.u[i] < (ll1 - ll0)){
      if(whichParam[i] == 1){
        s0 = s1
      }else{
        pi0 = pi1
      }
      ll0 = ll1
    }

    if(i > nBurnIn){
      chain[i - nBurnIn, ] = c(pi0, s0)
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
  fit$par = apply(chain, 2, mean)


  fitted = d.one.inflated.zeta(nvals, shape = fit$par[2], p = fit$par[1])
  names(fitted) = if(x$type == 'P'){
    paste0("P", nvals - 1)
  }else{
    paste0("S", nvals)
  }

  result = list(
    psData = x,
    fit = fit,
    pi = fit$par[1],
    shape =  fit$par[2],
    var.cov = cov(chain),
    fitted = fitted,
    chain = chain,
    model = "ziz",
    method = "bayes"
  )


  class(result) = "psFit"

  return(result)
}
