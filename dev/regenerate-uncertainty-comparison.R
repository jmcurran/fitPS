#!/usr/bin/env Rscript

# Regenerate the expensive uncertainty-comparison vignette assets.
#
# This script deliberately lives under dev/ so normal package builds do not
# rerun the 1,000-replicate ordinary and Rubin Bayesian bootstraps. It also
# verifies that plotUncertainty() reproduces an unconstrained Hscv/KDE contour
# construction from exactly the same ordinary-bootstrap replicates. This guards
# against reintroducing the positive-support transformation that distorted the
# joint bootstrap density geometry.

options(warn = 2)

devtools::load_all(".", quiet = TRUE)

data("Psurveys", package = "fitPS")
roux = Psurveys$roux
uncertaintyLevels = c(0.80, 0.95)

mleFit = fit(
  roux,
  model = zizModel(),
  nterms = 6
)

bayesFit = fit(
  roux,
  model = zizModel(),
  nterms = 6,
  method = "bayes",
  bayesOptions = list(posteriorMethod = "numerical")
)

bootFit = fit(
  roux,
  model = zizModel(),
  method = "bootstrap",
  nterms = 6,
  B = 1000,
  seed = 1234,
  silent = TRUE,
  parallel = FALSE
)

bayesBoot = fit(
  roux,
  model = zizModel(),
  method = "bayesianBootstrap",
  nterms = 6,
  B = 1000,
  seed = 1234
)

collectCoordinate = function(x, coordinate) {
  if (is.data.frame(x) && coordinate %in% names(x)) {
    return(x[[coordinate]])
  }
  if (is.list(x) && !is.null(x[[coordinate]]) && is.numeric(x[[coordinate]])) {
    return(x[[coordinate]])
  }
  if (is.list(x)) {
    return(unlist(
      lapply(x, collectCoordinate, coordinate = coordinate),
      use.names = FALSE
    ))
  }
  numeric(0)
}

referenceBootstrapRegion = function(replicates, level) {
  values = as.matrix(replicates[, c("pi", "shape"), drop = FALSE])
  values = values[complete.cases(values), , drop = FALSE]
  bandwidth = ks::Hscv(values)
  densityEstimate = ks::kde(
    x = values,
    H = bandwidth
  )
  contourHeights = ks::contourLevels(
    densityEstimate,
    cont = sort(100 * level),
    approx = TRUE
  )
  contours = grDevices::contourLines(
    x = densityEstimate$eval.points[[1L]],
    y = densityEstimate$eval.points[[2L]],
    z = densityEstimate$estimate,
    levels = contourHeights
  )
  list(
    density = densityEstimate,
    heights = contourHeights,
    contours = contours
  )
}

canonicalContour = function(contour) {
  points = cbind(x = contour$x, y = contour$y)
  startIndex = which.min(points[, "x"] + points[, "y"] * 1e-8)
  indices = c(
    seq.int(startIndex, nrow(points)),
    if (startIndex > 1L) seq_len(startIndex - 1L) else integer(0)
  )
  points[indices, , drop = FALSE]
}

compareContourSets = function(modernContours, referenceContours, tolerance = 1e-10) {
  if (length(modernContours) != length(referenceContours)) {
    stop("Modern and reference unconstrained KDE constructions produced different numbers of contours")
  }

  modern = lapply(modernContours, function(contour) {
    list(
      level = contour$densityLevel,
      points = canonicalContour(list(x = contour$x, y = contour$y))
    )
  })
  reference = lapply(referenceContours, function(contour) {
    list(
      level = contour$level,
      points = canonicalContour(contour)
    )
  })

  modernOrder = order(vapply(modern, `[[`, numeric(1L), "level"))
  referenceOrder = order(vapply(reference, `[[`, numeric(1L), "level"))
  modern = modern[modernOrder]
  reference = reference[referenceOrder]

  for (index in seq_along(modern)) {
    if (!isTRUE(all.equal(modern[[index]]$level, reference[[index]]$level, tolerance = tolerance))) {
      stop("Modern and reference unconstrained KDE contour heights differ")
    }
    if (!identical(dim(modern[[index]]$points), dim(reference[[index]]$points))) {
      stop("Modern and reference unconstrained KDE contours contain different numbers of coordinates")
    }

    forwardDifference = max(abs(modern[[index]]$points - reference[[index]]$points))
    reverseDifference = max(abs(
      modern[[index]]$points - reference[[index]]$points[nrow(reference[[index]]$points):1L, , drop = FALSE]
    ))
    if (min(forwardDifference, reverseDifference) > tolerance) {
      stop("plotUncertainty() does not preserve the reference unconstrained KDE contour geometry")
    }
  }

  invisible(TRUE)
}

expandRange = function(x, fraction = 0.05) {
  limits = range(x, finite = TRUE)
  width = diff(limits)
  if (!is.finite(width) || width <= 0) {
    width = max(abs(limits), 1)
  }
  limits + c(-1, 1) * fraction * width
}

# Capture the modern ordinary-bootstrap geometry without displaying it.
temporaryDevice = tempfile(fileext = ".pdf")
grDevices::pdf(temporaryDevice)
modernBootstrap = plotUncertainty(
  bootFit,
  level = uncertaintyLevels,
  main = "Ordinary bootstrap"
)
grDevices::dev.off()
unlink(temporaryDevice)

referenceBootstrap = referenceBootstrapRegion(
  bootFit$bootstrap$replicates,
  level = uncertaintyLevels
)
compareContourSets(
  modernContours = modernBootstrap$contours,
  referenceContours = referenceBootstrap$contours
)

principalAxisAngle = function(x, y) {
  covariance = stats::cov(cbind(x, y))
  vectors = eigen(covariance, symmetric = TRUE)$vectors
  (atan2(vectors[2L, 1L], vectors[1L, 1L]) * 180 / pi) %% 180
}

axisAngleDifference = function(first, second) {
  difference = abs(first - second) %% 180
  min(difference, 180 - difference)
}

bootstrapValues = as.data.frame(bootFit$bootstrap$replicates)
bootstrapValues = bootstrapValues[
  complete.cases(bootstrapValues[, c("pi", "shape")]),
  c("pi", "shape"),
  drop = FALSE
]
cloudAngle = principalAxisAngle(bootstrapValues$pi, bootstrapValues$shape)
outerContours = Filter(
  function(contour) isTRUE(all.equal(contour$level, max(uncertaintyLevels))),
  modernBootstrap$contours
)
if (length(outerContours) < 1L) {
  stop("The regenerated ordinary-bootstrap plot did not contain the requested outer contour")
}
outerSizes = vapply(outerContours, function(contour) length(contour$x), integer(1L))
outerContour = outerContours[[which.max(outerSizes)]]
contourAngle = principalAxisAngle(outerContour$x, outerContour$y)
orientationDifference = axisAngleDifference(cloudAngle, contourAngle)
if (orientationDifference > 3.5) {
  stop(sprintf(
    "Ordinary-bootstrap 95%% KDE contour is %.2f degrees away from the replicate-cloud principal axis",
    orientationDifference
  ))
}

mleRegions = confint(
  mleFit,
  level = uncertaintyLevels,
  silent = TRUE
)
bayesRegions = credint(
  bayesFit,
  level = uncertaintyLevels,
  silent = TRUE
)
bootstrapReplicates = bootFit$bootstrap$replicates
bayesianBootstrapReplicates = bayesBoot$replicates$parameters
bootstrapContours = referenceBootstrap$contours
bayesianBootstrapRegion = referenceBootstrapRegion(
  bayesianBootstrapReplicates,
  level = uncertaintyLevels
)
bayesianBootstrapContours = bayesianBootstrapRegion$contours

sharedPi = c(
  collectCoordinate(mleRegions, "pi"),
  bootstrapReplicates$pi,
  collectCoordinate(bootstrapContours, "x"),
  collectCoordinate(bayesRegions, "pi"),
  bayesianBootstrapReplicates$pi,
  collectCoordinate(bayesianBootstrapContours, "x")
)
sharedShape = c(
  collectCoordinate(mleRegions, "shape"),
  bootstrapReplicates$shape,
  collectCoordinate(bootstrapContours, "y"),
  collectCoordinate(bayesRegions, "shape"),
  bayesianBootstrapReplicates$shape,
  collectCoordinate(bayesianBootstrapContours, "y")
)
sharedXlim = expandRange(sharedPi)
sharedYlim = expandRange(sharedShape)

figuresDirectory = file.path("vignettes", "figures")
dir.create(figuresDirectory, recursive = TRUE, showWarnings = FALSE)
figurePath = file.path(figuresDirectory, "uncertainty-comparison.png")
containmentPath = file.path(figuresDirectory, "uncertainty-containment.csv")
cachePath = file.path(figuresDirectory, "uncertainty-bootstrap-cache.rds")

png(
  filename = figurePath,
  width = 1800,
  height = 1800,
  res = 200
)
oldPar = par(no.readonly = TRUE)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

plotUncertainty(
  mleFit,
  level = uncertaintyLevels,
  xlim = sharedXlim,
  ylim = sharedYlim,
  main = "Profile likelihood"
)
bootstrapPlotInfo = plotUncertainty(
  bootFit,
  level = uncertaintyLevels,
  xlim = sharedXlim,
  ylim = sharedYlim,
  main = "Ordinary bootstrap"
)
plotUncertainty(
  bayesFit,
  level = uncertaintyLevels,
  xlim = sharedXlim,
  ylim = sharedYlim,
  main = "Parametric Bayes"
)
bayesianBootstrapPlotInfo = plotUncertainty(
  bayesBoot,
  level = uncertaintyLevels,
  xlim = sharedXlim,
  ylim = sharedYlim,
  main = "Bayesian Bootstrap"
)

par(oldPar)
dev.off()

containmentTable = rbind(
  data.frame(
    method = "Ordinary bootstrap",
    nominal = bootstrapPlotInfo$containment$level,
    empirical = bootstrapPlotInfo$containment$containment
  ),
  data.frame(
    method = "Bayesian Bootstrap",
    nominal = bayesianBootstrapPlotInfo$containment$level,
    empirical = bayesianBootstrapPlotInfo$containment$containment
  )
)
utils::write.csv(containmentTable, containmentPath, row.names = FALSE)
saveRDS(
  list(
    bootFit = bootFit,
    bayesBoot = bayesBoot
  ),
  cachePath
)

cat(sprintf(
  "Verified unconstrained ordinary-bootstrap KDE geometry; 95%% contour/cloud axis difference: %.2f degrees.\n",
  orientationDifference
))
cat("Wrote:", figurePath, "\n")
cat("Wrote:", containmentPath, "\n")
cat("Wrote:", cachePath, "\n")
