fitZIDist.old = function(x, nterms = 10,
                   start = c(0.5, 1),
                   ...){
  nvals = 1:nterms
  if(!is(x, "psData")){
    stop("x must be an object of class psData")
  }

  if(start[1] <= 0 || start[1] >= 1){
    stop("The starting value for pi must be in (0, 1)")
  }

  if(start[2] <= 0){
    stop("The Zeta function is undefined for shape = 0. Choose a start value > 0.")
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

    -sum(d.one.inflated.zeta(y, shape, p, log = TRUE))
  }

  # fit = nlminb(start = start,
  #              objective = logLik,
  #              lower = 1)

  fit = optim(par = start,
              fn = logLik,
              method = "L-BFGS-B",
              lower = c(0, 0),
              upper  = c(1 - .Machine$double.eps, Inf),
              hessian = TRUE)

  fitted = d.one.inflated.zeta(nvals, shape = fit$par[2], p = fit$par[1])
  names(fitted) = if(x$type == 'P'){
    paste0("P", nvals - 1)
  }else{
    paste0("S", nvals)
  }

  result = list(
    psData = x,
    fit = fit,
    shape =  fit$par[2],
    pi = fit$par[1],
    var.cov = solve(fit$hessian),
    fitted = fitted,
    zeroInflated = TRUE
  )


  class(result) = "psFit"

  return(result)
}

fitZIDist = function(x, nterms = 10,
                     start = c(0.5, 1),
                     ...){
  nvals = 1:nterms
  if(!is(x, "psData")){
    stop("x must be an object of class psData")
  }

  if(start[1] <= 0 || start[1] >= 1){
    stop("The starting value for pi must be in (0, 1)")
  }

  if(start[2] <= 0){
    stop("The Zeta function is undefined for shape = 0. Choose a start value > 0.")
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

    -sum(x$data$rn * log(rval))
  }

  # fit = nlminb(start = start,
  #              objective = logLik,
  #              lower = 1)

  fit = optim(par = start,
              fn = logLik,
              method = "L-BFGS-B",
              lower = c(0, 0),
              upper  = c(1 - .Machine$double.eps, Inf),
              hessian = TRUE)

  fitted = d.one.inflated.zeta(nvals, shape = fit$par[2], p = fit$par[1])
  names(fitted) = if(x$type == 'P'){
    paste0("P", nvals - 1)
  }else{
    paste0("S", nvals)
  }

  result = list(
    psData = x,
    fit = fit,
    shape =  fit$par[2],
    pi = fit$par[1],
    var.cov = solve(fit$hessian),
    fitted = fitted,
    zeroInflated = TRUE
  )


  class(result) = "psFit"

  return(result)
}

