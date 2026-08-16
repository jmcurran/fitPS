#!/usr/bin/env Rscript

# Regenerate the expensive uncertainty-comparison vignette assets.
#
# This script deliberately lives under dev/ so normal package builds do not
# rerun the 1,000-replicate ordinary and Rubin Bayesian bootstraps. It also
# verifies that plotUncertainty() reproduces the historical bootCI() KDE
# contour construction from exactly the same ordinary-bootstrap replicates.

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

legacyBootstrapRegion = function(replicates, level) {
  values = as.matrix(replicates[, c("pi", "shape"), drop = FALSE])
  values = values[complete.cases(values), , drop = FALSE]
  bandwidth = ks::Hscv(values)
  densityEstimate = ks::kde(
    x = values,
    H = bandwidth,
    positive = TRUE
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

compareContourSets = function(modernContours, legacyContours, tolerance = 1e-10) {
  if (length(modernContours) != length(legacyContours)) {
    stop("Modern and historical KDE constructions produced different numbers of contours")
  }

  modern = lapply(modernContours, function(contour) {
    list(
      level = contour$densityLevel,
      points = canonicalContour(list(x = contour$x, y = contour$y))
    )
  })
  legacy = lapply(legacyContours, function(contour) {
    list(
      level = contour$level,
      points = canonicalContour(contour)
    )
  })

  modernOrder = order(vapply(modern, `[[`, numeric(1L), "level"))
  legacyOrder = order(vapply(legacy, `[[`, numeric(1L), "level"))
  modern = modern[modernOrder]
  legacy = legacy[legacyOrder]

  for (index in seq_along(modern)) {
    if (!isTRUE(all.equal(modern[[index]]$level, legacy[[index]]$level, tolerance = tolerance))) {
      stop("Modern and historical KDE contour heights differ")
    }
    if (!identical(dim(modern[[index]]$points), dim(legacy[[index]]$points))) {
      stop("Modern and historical KDE contours contain different numbers of coordinates")
    }

    forwardDifference = max(abs(modern[[index]]$points - legacy[[index]]$points))
    reverseDifference = max(abs(
      modern[[index]]$points - legacy[[index]]$points[nrow(legacy[[index]]$points):1L, , drop = FALSE]
    ))
    if (min(forwardDifference, reverseDifference) > tolerance) {
      stop("plotUncertainty() does not preserve the historical bootCI() KDE contour geometry")
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

legacyBootstrap = legacyBootstrapRegion(
  bootFit$bootstrap$replicates,
  level = uncertaintyLevels
)
compareContourSets(
  modernContours = modernBootstrap$contours,
  legacyContours = legacyBootstrap$contours
)

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
bootstrapContours = legacyBootstrap$contours
bayesianBootstrapRegion = legacyBootstrapRegion(
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

cat("Verified plotUncertainty() ordinary-bootstrap contours against the historical bootCI() KDE construction.\n")
cat("Wrote:", figurePath, "\n")
cat("Wrote:", containmentPath, "\n")
cat("Wrote:", cachePath, "\n")
