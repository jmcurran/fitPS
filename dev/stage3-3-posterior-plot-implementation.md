# Stage 3.3 posterior plotting implementation

## Purpose

Stage 3.3 implements a first posterior-density plotting API for Bayesian fitPS model fits.

## User-facing change

The stage adds `plotPosterior()`, an exported plotting helper for `psFit` objects. It keeps the existing `plot.psFit()` behavior unchanged, so fitted probability plots remain backward compatible.

## Supported fit objects

The implementation supports:

- zeta Bayesian MCMC fits with a numeric `chain`;
- zeta numerical-integration fits with a stored `pdf` function;
- zero-inflated Bayesian MCMC fits with a data-frame or matrix `chain` containing `shape` and `pi` columns.

The default parameter is `shape`. Zero-inflated Bayesian fits can also use `parameter = "pi"`.

## Plot behavior

The plot shows the marginal posterior density, the stored posterior point estimate, and an equal-tail credible interval by default. For MCMC fits, the density and interval are computed from stored samples. For numerical-integration fits, the stored posterior density is evaluated on an automatically chosen grid, converted to a linear interpolation function with `approxfun()`, normalized with `integrate()`, and inverted with `uniroot()` for equal-tail interval construction. This avoids bespoke trapezoid helpers and keeps the numerical calculation faithful to the stored grid.

## Testing approach

Tests use small fake `psFit` objects so they remain deterministic, offline, and fast. They cover zeta MCMC samples, zero-inflated `shape` and `pi` samples, numerical posterior-density functions, and error handling for unsupported fits.
