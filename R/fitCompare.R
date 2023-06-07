fitCompare = function(x, start = list(zeta = 1, ziz = c(0.5, 1)), ...){
  fit.zeta = fitDist(x, start = start$zeta, ...)
  fit.ziz = fitZIdist(x, start = start$ziz, ...)
  p.zeta = probfun(fit.zeta)
  p.ziz = probfun(fit.ziz)
  raw = x$data$rn / sum(x$data$rn)

  nmax = max(x$data$n) + 1

  if(x$type == "P"){
    fitted = data.frame(raw = rep(0, nmax + 1),
                        zeta = rep(0, nmax + 1),
                        ziz = rep(0, nmax + 1))
    fitted$raw[x$data$n + 1] = raw
    fitted$zeta = p.zeta(0:nmax)
    fitted$ziz = p.ziz(0:nmax)
  }else{
    fitted = data.frame(raw = rep(0, nmax),
                        zeta = rep(0, nmax),
                        ziz = rep(0, nmax))
    fitted$raw[x$data$n] = raw
    fitted$zeta = p.zeta(1:nmax)
    fitted$ziz = p.ziz(1:nmax)
  }

  print(fitted)
  invisible(list(raw = raw, fit.zeta = fit.zeta, fit.ziz = fit.ziz, fitted = fitted))
}
