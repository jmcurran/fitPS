#' Bootstrap confidence intervals or regions
#'
#' Use boostrapping to generate confidence intervals, or
#' confidence regions in the case of the zero-inflated model
#'
#' @param x a object of class \code{psFit}---see \code{\link{readDta}} for more
#' details.
#' @param B the number of bootstrap samples to take.
#' @param model which model to fit to the data, either \code{"zeta"} or
#' \code{"zi.zeta}. Maybe abbreviated to \code{"z"} and \code{"zi"}. Default is
#' \code{"zeta"}.
#' @param level the confidence level required---restricted to [0.75, 1)
#'
#' importFrom doParallel regusterDoParallel
#' import foreach
#' importFrom pbapply pblapply pbsapply pboptions
#' importFrom parallel detectCores makeCluster parApply parLapply parSapply stopCluster
bootCI = function(x,
                  level = 0.95,
                  B = 2000,
                  model = c("zeta", "zi.zeta"),
                  silent = FALSE,
                  parallel = TRUE,
                  progressBar = FALSE,
                  pbopts = list(type = "txt"){

}


bootFit = function(x, B = 2000, model = c("zeta", "zi.zeta"),
                   silent = FALSE,
                   parallel = TRUE,
                   progressBar = FALSE,
                   pbopts = list(type = "txt")){
  yvals = rep(x$data$n, x$data$rn)
  n = length(yvals)

  to.psData = function(y, type){
    tbl = table(y)
    counts = as.vector(tbl)

    r = list(data = data.frame(n =  as.numeric(names(tbl)),
                               rn =  counts))
    r$type = type
    class(r) = "psData"

    return(r)
  }

  environment(to.psData) = baseenv()

  boot.y = matrix(sample(yvals, n * B, replace = TRUE), nrow = B)

  model = match.arg(model)

  ## just evaluate the progress bar options once instead of calling it everywhere
  if(progressBar){
    opb = do.call(get("pboptions", asNamespace("pbapply")), pbopts)
    on.exit(pbapply::pboptions(opb))
  }

  if(parallel){
    ncores = parallel::detectCores()
    cl = parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    if(!silent){
      cat("Creating bootstrapped data sets\n")
    }

    boot.y = if(progressBar){
      pbapply::pbapply(X = boot.y, MARGIN = 1, FUN = to.psData, type = x$type, cl = cl)
    }else{
      ##parallel::parApply(cl = cl, X = boot.y, MARGIN = 1, FUN = to.psData, type = x$type)
      foreach(i = 1:nrow(boot.y)) %dopar% to.psData(boot.y[i,], type = x$type)
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
        # parallel::parSapply(cl = cl, X = boot.y, FUN = function(y){
        #   fitDist(y)$shape
        # })
        foreach(i = seq_along(boot.y), .combine = 'c') %dopar% {
          r = fitDist(boot.y[[i]])
          r$shape
        }
      }
    }else{
      results = if(progressBar){
        pbapply::pblapply(X = boot.y, fun  = function(y){
          fit = fitZIDist(y)
          return(c(fit$pi, fit$shape))
        }, cl = cl)
      }else{
        # parallel::parLapply(cl = cl, X = boot.y, fun  = function(y){
        #   fit = fitZIDist(y)
        #   return(c(fit$pi, fit$shape))
        # })
        foreach(i = seq_along(boot.y)) %dopar% {
          fit = fitZIDist(boot.y[[i]])
          c(fit$pi, fit$shape)
        }
      }
      results = as.data.frame(do.call("rbind", results))
      names(results) = c("pi", "shape")
    }
    parallel::stopCluster(cl)
  }else{
    boot.y = if(progressBar){
        pbapply::pbapply(boot.y, 1, to.psData, type = x$type)
      }else{
        apply(boot.y, 1, to.psData, type = x$type)
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
      results = as.data.frame(do.call("rbind", results))
      names(results) = c("pi", "shape")
    }
  }

  return(results)
}