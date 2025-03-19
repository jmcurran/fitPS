zetaloglikelihood = function(x, shape){
  offset = ifelse(x$type == 'P', 1, 0)
  if(length(shape) > 1){
    ll = function(s){
      sum(x$data$rn * VGAM::dzeta(x$data$n + offset, shape = s, log = TRUE))
    }
    sapply(shape, ll)
  }else{
    sum(x$data$rn * VGAM::dzeta(x$data$n + offset, shape = shape, log = TRUE))
  }
}

zll = zetaloglikelihood

zetalikelihood = function(x, shape){
  exp(zll(x, shape))
}

zetaunifprior = function(s, a = exp(-2), b = exp(2)){
  result = numeric(length(s))
  inRange = s >= a & s <= b
  result[inRange] = 1 / ((b - a) * s[inRange])
  result[!inRange] = 0
  return(result)
}

zup = zetaunifprior
