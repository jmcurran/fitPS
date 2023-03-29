plZIZ  = function(x, level = 0.95,
                 grid.Pi = seq(0.5, 1 - .Machine$double.eps, length = 100),
                 grid.Shape = seq(1, x$shape + 4 * sx, by = 0.01),
                 silent = FALSE){
  if(!x$zeroInflated){
    stop("This function is for the ZIZ model only")
  }

  if(any(level <= 0 | level >= 1)){
    stop("The elements of level must be in the interval (0, 1)")
  }

  sx = sqrt(diag(x$var.cov))[2]
  l0 = -x$fit$val
  dataf = x$psData$data
  if(x$psData$type == "P"){
    dataf$n = dataf$n + 1
  }

  if(!silent){
    message("Computing contours. This may take a few seconds.")
  }

  r = lapply(grid.Pi, function(p){
    sapply(grid.Shape, function(si){
      -2 * (zi.loglik(dataf, c(p, si)) - l0)
    })
  })

  r = do.call("rbind", r)

  cr = lapply(level, function(l){
    qstar2 = qchisq(l, 2)
    r1 = r - qstar2
    confRegion = contourLines(grid.Pi, grid.Shape, r1, levels = 0)[[1]]
    confRegion = data.frame(pi = confRegion$x, shape = confRegion$y)
  })

  names(cr) = as.character(level)

  return(cr)
}
