#' importFrom pbapply pblapply pbsapply pboptions
#' importFrom parallel detectCores makeCluster parApply, parLapply parSapply stopCluster
bootFit = function(x, B = 2000, model = c("zeta", "zi.zeta"),
                   silent = FALSE,
                   parallel = TRUE,
                   progressBar = FALSE,
                   pbopts = list(type = "txt")){
  yvals = rep(x$data$n, x$data$rn)
  n = length(yvals)

  to.psData = function(y){
    tbl = table(y)
    counts = as.vector(tbl)

    r = x
    r$data = data.frame(n =  as.numeric(names(tbl)),
                        rn =  counts)

    return(r)
  }

  boot.y = matrix(sample(yvals, n * B, replace = TRUE), nrow = B)

  model = match.arg(model)

  ## just evaluate the progress bar options once instead of calling it everywhere
  if(progressBar){
    obp = do.call(get("pboptions", asNamespace("pbapply")), pbopts)
    on.exit(pbapply::pboptions(opb))
  }

  if(parallel){
    ncores = parallel::detectCores()
    cl = parallel::makeCluster(ncores)

    if(!silent){
      cat("Creating bootstrapped data sets\n")
    }

    boot.y = if(progressBar){
      pbapply::pbapply(X = boot.y, MARGIN = 1, FUN = to.psData, cl = cl)
    }else{
      parallel::parApply(cl = cl, X = boot.y, MARGIN = 1, FUN = to.psData)
    }

    if(!silent){
      cat("Estimating parameters for each bootstrapped data set\n")
    }

    if(model == "zeta"){
      results = if(progressBar){
        pbapply::pbsapply(X = boot.y, FUN = function(y){
          fitDist(y)$shape
        }, cl = cl)
      }else{
        parallel::parSapply(cl = cl, X = boot.y, FUN = function(y){
          fitDist(y)$shape
        })
      }
    }else{
      results = if(progressBar){
        pbapply::pblapply(X = boot.y, fun  = function(y){
          fit = fitZIDist(y)
          return(c(fit$pi, fit$shape))
        }, cl = cl)
      }else{
        parallel::parLapply(cl = cl, X = boot.y, fun  = function(y){
          fit = fitZIDist(y)
          return(c(fit$pi, fit$shape))
        })
      }
      results = do.call("rbind", results)
      names(results) = c("pi", "shape")
    }
    parallel::stopCluster(cl)
  }else{
    boot.y = if(progressBar){
        pbapply::apply(boot.y, 1, to.psData)
      }else{
        apply(boot.y, 1, to.psData)
      }
    if(model == "zeta"){
      results = if(progressBar){
        pbapply::pbsapply(boot.y, function(y){
          fitDist(y)$shape
        })
      }else{
        sapply(boot.y, function(y){
          fitDist(y)$shape
        })
      }
    }else{
      results = if(progressBar){
        pbapply::pblapply(boot.y, function(y){
          fit = fitZIDist(y)
          return(c(fit$pi, fit$shape))
        })
      }else{
        lapply(boot.y, function(y){
          fit = fitZIDist(y)
          return(c(fit$pi, fit$shape))
        })
      }
      results = do.call("rbind", results)
      names(results) = c("pi", "shape")
    }
  }

  return(results)
}
