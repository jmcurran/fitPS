# Stage 9.6.3 mode-centred cubature repair

Stage 9.6.2 integrated directly over finite natural-parameter bounds. For the
external Poisson-normal proof model, adaptive cubature could still miss a
narrow posterior region and return a zero normalization integral.

Stage 9.6.3 keeps the same finite model-supplied natural bounds, but maps each
parameter interval piecewise to the unit interval with the optimized posterior
mode fixed at 0.5. The two-dimensional posterior mode is therefore located at
(0.5, 0.5), while the piecewise-linear Jacobian preserves the exact integral
over the original natural-parameter rectangle.

This is a numerical-stability repair only. The Stage 9.6 policy remains:

- one parameter: base `integrate()`;
- two parameters: adaptive `cubature::hcubature()`;
- three or more parameters: use MCMC.
