# Stage 9.5: generic one-dimensional numerical posterior engine

Stage 9.5 makes the existing numerical posterior strategy genuinely generic for
one-parameter `psModel` objects.

## Architecture

The model owns the likelihood and prior mathematics through
`modelLogLikelihood()` and `modelLogPrior()`. For the numerical engine,
`modelBayesControl()` also supplies named natural-scale `start`, `lower`, and
`upper` values. Bounds may be infinite, so external models are not restricted
to the legacy finite `psPrior` range.

The numerical engine now owns posterior optimization, normalization, posterior
moments, the stored scalar density function, and generic derived-probability
summaries. This keeps the model/engine boundary consistent with the generic
MCMC path while preserving natural-scale density functions used by plotting and
DIC.

## Built-in models

The zeta and logarithmic model-specific numerical fitter and numerical
probability-summary methods are removed. Both models use the shared generic
one-dimensional engine and retain their model-specific likelihood, prior,
parameter validation, and Bayesian controls.

ZIZ remains deliberately model-specific for numerical fitting because its
established posterior is a two-dimensional grid. Stage 9.5 does not alter its
numerical, Laplace, or importance implementations.

## External proof

The external Poisson test model now declares both numerical and MCMC support.
Its Gamma prior uses natural-scale bounds `(0, Inf)`, demonstrating that the
public numerical contract is not tied to `psPrior`. A deterministic regression
test compares the generic numerical posterior mean and variance with the
conjugate Gamma-Poisson posterior and checks density normalization, fitted
probabilities, and DIC.
