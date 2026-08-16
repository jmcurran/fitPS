# Stage 9.6.1 test parse repair

Stage 9.6 failed before the test suite could exercise the new two-dimensional numerical posterior engine because `tests/testthat/test-generic-numerical-posterior.R` ended its final `test_that()` call with `}` instead of `})`.

This repair changes only that test syntax. The Stage 9.6 implementation remains unchanged: one-dimensional numerical posterior fitting uses `integrate()`, two-dimensional fitting uses `cubature::hcubature()`, and models with more than two parameters are directed to MCMC.

Build `1.1.1.012` was consumed by the failed Stage 9.6 attempt. Stage 9.6.1 therefore advances to build `1.1.1.013`.
