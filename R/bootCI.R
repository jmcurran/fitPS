#' Bootstrap confidence intervals or regions
#'
#' Use boostrapping to generate confidence intervals, or confidence regions in
#' the case of the zero-inflated model.
#'
#' @param x a object of class \code{psFit}---see \code{\link{readData}} for more
#'   details.
#' @param level the confidence level required---restricted to [0.75, 1). This
#'   may be a vector, in which case multiple intervals, or confidence regions
#'   will be returned.
#' @param B the number of bootstrap samples to take.
#' @param model which model to fit to the data, either \code{"zeta"} or
#'   \code{"zi.zeta"}. Maybe abbreviated to \code{"z"} and \code{"zi"}. Default
#'   is \code{"zeta"}.
#' @param silent if \code{TRUE}, then no output will be displayed whilst the
#'   bootstrapping is being undertaken. \code{plot} if \code{TRUE} then the
#'   contours for the confidence region will be plotted. This only works if
#'   \code{model = "zi.zeta"}. It is ignored otherwise. \code{parallel} if
#'   \code{TRUE} then the bootstrapping is performed in parallel.
#' @param plot if \code{TRUE} and \code{model == "zi.zeta"}, then a plot of the
#'   bootstrapped values will be produced and confidence contour lines will be
#'   drawn for each value in level.
#' @param parallel if \code{TRUE}, then the package will attempt to use multiple
#'   cores to speed up computation.
#' @param progressBar if \code{TRUE}, then progress bars will be displayed to
#'   show progress on the bootstrapping.
#' @param pbopts a list of arguments for the \code{\link[pbapply]{pboptions}}
#'   function that affect the progress bars. Ignored if \code{progressBar =
#'   FALSE}.
#' @param ... other arguments.
#'
#' @details This function uses bootstrapping to compute a confidence interval
#'   for the shape parameter in the case of the zeta model and a confidence
#'   region in the case of the zero-inflated zeta model. A smoothed bootstrap
#'   approach is taken rather than a simple percentile method. The kernel
#'   density estimation is performed by the \code{ks} package using a smoothed
#'   cross-validated bandwidth selection procedure.
#'
#' @returns If \code{model == "zeta"}, then either a \code{vector} or a
#'   \code{data.frame} with elements/columns named \code{"lower"} and
#'   \code{"upper"} representing the lower and upper bounds of the confidence
#'   interval(s). Multiple bounds are returned in a \code{data.frame} when
#'   \code{level} has more than one value. If \code{model == "zi.zeta"}, then a
#'   list with length equal to the length of \code{level} is returned. The name
#'   of each element in the list is the level with % attached. For example if
#'   \code{level == 0.95}, then the list has a single element named
#'   \code{"95\%"}. Each element of the list consists of a \code{list} with
#'   elements named \code{pi} and \code{shape} which specify the coordinates of
#'   the contour for that level. There is a third element named \code{level}
#'   which gives the height of the kernel density estimate at that contour.
#'
#' @examples
#' \dontrun{
#' data(Psurveys)
#' roux = Psurveys$roux
#' confRegion = bootCI(roux, model = "zi.zeta", parallel = FALSE, plot = TRUE)
#'
#' ## This will not work unless you have the sp package installed
#' ## Count how many of the points lie within the 95% confidence region
#' table(sp::point.in.polygon(fit$pi,fit$shape, confRegion$pi, confRegion$shape))
#' }
#' @importFrom doParallel registerDoParallel
#' @import foreach
#' @importFrom grDevices contourLines
#' @importFrom graphics polygon
#' @importFrom iterators iter
#' @importFrom ks contourLevels kcde kde Hscv
#' @importFrom pbapply pblapply pbsapply pboptions
#' @importFrom parallel detectCores makeCluster parApply parLapply parSapply
#'   stopCluster
#' @importFrom stats approxfun
#' @export
bootCI = function(x, ...){
  UseMethod("bootCI", x)
}

#' @describeIn bootCI Bootstrap confidence intervals or regions
#' @export
bootCI.default = function(x,
                          level = 0.95,
                          B = 2000,
                          model = c("zeta", "zi.zeta"),
                          silent = FALSE,
                          plot = FALSE,
                          parallel = TRUE,
                          progressBar = FALSE,
                          pbopts = list(type = "txt"),
                          ...){

  model = match.arg(model)

  if(any(level < 0.75 | level >= 1)){
    stop("The entries level must be values between 0.75 and 1 (not inclusive).")
  }

  fit = bootFit(x = x,
                B = B,
                model = model,
                silent = silent,
                parallel = parallel,
                progressBar = progressBar,
                pbopts = pbopts)


  if(model == "zeta"){
    ## estimate bandwidth
    h = ks::hscv(fit)

    fhat = ks::kcde(fit, h) ## computes the CDF based on the KDE
    FxInv = approxfun(fhat$estimate, fhat$eval.points)

    if(length(level) == 1){
      alpha2 = 0.5 * (1 - level)
      ci = c(FxInv(alpha2), FxInv(1 - alpha2))
      names(ci) = c("lower", "upper")
    }else{
      level = sort(level)
      alpha2 = 0.5 * (1 - level)
      ci = data.frame(lower = FxInv(alpha2),
                      upper = FxInv(1 - alpha2))
    }
    return(ci)
  }else{
    ## estimate bandwidth
    H = ks::Hscv(fit)

    fhat = ks::kde(fit, H)
    cont = sort(100 * level)
    confRegion = contourLines(x = fhat$eval.points[[1]],
                 y = fhat$eval.points[[2]],
                 z = fhat$estimate,
                 levels = ks::contourLevels(fhat, cont = cont, approx = TRUE))
    names(confRegion) = paste0(cont,"%")
    confRegion = lapply(confRegion, function(l){
      names(l)[2:3] = c("pi", "shape")
      return(l)
    })

    if(plot){
      plot(fit, pch = 'x', col = 'grey')
      for(cr in seq_along(confRegion)){
        polygon(confRegion[[cr]]$pi, confRegion[[cr]]$shape, border = "red", lwd = 2)
      }
    }

    #names(confRegion) = c("pi", "shape")
    return(confRegion)
  }
}

#' @describeIn bootCI Bootstrap confidence intervals or regions
#' @export
bootCI.psData = function(x, ...){
  return(bootCI.default(x = x, ...))
}

#' @describeIn bootCI Bootstrap confidence intervals or regions
#' @export
bootCI.psFit = function(x, model = ifelse(!x$zeroInflated, "zeta", "zi.zeta"),
                        ...){
  return(bootCI.default(x = x$psData, model = model, ...))
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
      foreach(row = iterators::iter(boot.y, by = "row")) %dopar% {to.psData(row, type = x$type)}
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
        foreach(x = boot.y, .combine = 'c') %dopar% {
          r = fitDist(x)
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
        foreach(x = boot.y) %dopar% {
          fit = fitZIDist(x)
          c(fit$pi, fit$shape)
        }
      }
      results = as.data.frame(do.call("rbind", results))
      names(results) = c("pi", "shape")
    }
    parallel::stopCluster(cl)
  }else{
    if(!silent){
      cat("Creating bootstrapped data sets\n")
    }

    boot.y = if(progressBar){
      pbapply::pbapply(boot.y, 1, to.psData, type = x$type)
    }else{
      apply(boot.y, 1, to.psData, type = x$type)
    }

    if(!silent){
      cat("Estimating parameters for each bootstrapped data set\n")
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

