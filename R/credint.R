#' Bayesian credible intervals or regions
#'
#' Use kernel density estimation to generate credible intervals, or credible
#' regions in the case of the zero-inflated model.
#'
#' @param psFit a object of class \code{psFit}.
#' @param level the credible level required---restricted to [0.75, 1). This may
#'  be a vector, in which case multiple intervals, or credible regions will be
#'  returned.
#' @param plot if \code{TRUE} and \code{model == "ziz"}, then a plot of the
#'  bootstrapped values will be produced and confidence contour lines will be
#'  drawn for each value in level.
#' @param silent if \code{TRUE}, then no output will be displayed whilst the the
#'  kernel density estimation is being undertaken.
#' @param ... other arguments fed to plot. If \code{plot == FALSE}, then these
#'  will be ignored
#'
#' @aliases credInt
#' @details This function uses kernel density estimation to compute a Bayesian
#'  credible interval for the shape parameter in the case of the zeta model and
#'  a credible region in the case of the zero-inflated zeta model. A smoothing
#'  approach is taken rather than a simple percentile method. The
#'  kernel density estimation is performed by the \code{ks} package using a
#'  smoothed cross-validated bandwidth selection procedure.
#'
#' @returns
#'  If \code{psData$model == "zeta"}, then either a \code{vector} or a \code{data.frame}
#'  with elements/columns named \code{"lower"} and \code{"upper"} representing
#'  the lower and upper bounds of the confidence interval(s). Multiple bounds
#'  are returned in a \code{data.frame} when \code{level} has more than one
#'  value. If \code{psData$model == "ziz"}, then a list with length equal to the
#'  length of \code{level} is returned. The name of each element in the list is
#'  the level with % attached. For example if \code{level == 0.95}, then the
#'  list has a single element named \code{"95\%"}. It is possible for there to
#'  be multiple contours for the confidence region for a given \code{level}. If
#'  there is only one contour for each value of \code{level}, then each element
#'  of the list consists of a \code{list} with elements named \code{pi} and
#'  \code{shape} which specify the coordinates of the contour(s) for that level.
#'  There is a third element named \code{level} which gives the height of the
#'  kernel density estimate at that contour. If there are multiple contours for
#'  a given value of \code{level} then each list element is a list of lists with
#'  the structure given above (\code{level}, \code{pi}, and \code{shape}). NOTE:
#'  it is quite possible that there are multiple contours for a given height. If
#'  you want a way of thinking about this consider a mountain range with two
#'  mountains of equal height. If you draw the contours for (almost) any
#'  elevation, then you would expect to capture a region from each mountain.
#'
#' @examples
#' \dontrun{
#' data(Psurveys)
#' roux = Psurveys$roux
#' fit = fitzidist(roux, method == "bayes")
#' credRegion = credint(roux, plot = TRUE)
#'
#' ## This will not work unless you have the sp package installed
#' ## Count how many of the points lie within the 95% confidence region
#' lapply(credRegion, function(cr){
#'   table(sp::point.in.polygon(fit$pi,fit$shape, cr$pi, cr$shape))
#'. })
#' }
#' @export
credint = function(psFit,
                   level = 0.95,
                   plot = FALSE,
                   silent = FALSE,
                   ...){

  if(!is(psFit, "psFit")){
    stop("This function only works with objects of class psFit.\nYou must run fitDist or fitZIDist first.")
  }
  model = psFit$model
  if(psFit$method != "bayes"){
    stop("This method is for Bayesian estimation only.\nUse confint or bootCI instead.")
  }

  if(any(level < 0.75 | level >= 1)){
    stop("The entries level must be values between 0.75 and 1 (not inclusive).")
  }



  if(!silent){
    cat("Computing contours\n")
  }

  if(model == "zeta"){
    ## estimate bandwidth
    if(!silent){
      cat("\t-- Estimating bandwidth\n")
    }
    h = ks::hscv(psFit$chain)

    if(!silent){
      cat("\t-- Computing KCDE\n")
    }
    fhat = ks::kcde(psFit$chain, h) ## computes the CDF based on the KDE
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
    if(!silent){
      cat("\t-- Estimating bandwidth\n")
    }
    H = ks::Hscv(psFit$chain)

    if(!silent){
      cat("\t-- Computing KDE\n")
    }
    fhat = ks::kde(psFit$chain, H, positive = TRUE)
    cont = sort(100 * level)
    levels = ks::contourLevels(fhat, cont = cont, approx = TRUE)
    credRegion = contourLines(x = fhat$eval.points[[1]],
                 y = fhat$eval.points[[2]],
                 z = fhat$estimate,
                 levels = levels)

    credRegion = lapply(credRegion, function(l){
      names(l)[2:3] = c("pi", "shape")
      return(l)
    })

    ## This section reorganises credRegion into a list of lists, with each
    ## main element corresponding to a confidence level. This recognises that
    ## there may be more than one contour for each credible level

    cr.levels = sapply(credRegion, function(cr)cr$level)
    i = match(cr.levels, levels)
    credRegion = split(credRegion, i)
    names(credRegion) = paste0(cont,"%")

    if(plot){
      plot(psFit$chain, pch = 'x', col = 'grey', ...)
      for(l in seq_along(credRegion)){
        for(cr in seq_along(credRegion[[l]])){
          polygon(credRegion[[l]][[cr]]$pi, credRegion[[l]][[cr]]$shape, border = "red", lwd = 2)
        }
      }
    }

    ## If there is only one contour per level, then remove the list structure
    if(length(credRegion) > 1 && all(sapply(credRegion, length) == 1)){
      credRegion = lapply(credRegion, function(l)unlist(l, recursive = FALSE))
    }
    return(credRegion)
  }
}

#' @describeIn credint Bayesian credible intervals or regions
#' @export
credInt = credint
