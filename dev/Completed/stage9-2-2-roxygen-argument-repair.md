# Stage 9.2.2 roxygen argument repair

Stage 9.2.1 reached package checking but failed because roxygen-generated Rd files included method-specific arguments that were not documented on the shared help pages.

This repair documents `shape1` and `shape2` on the `modelLogPrior()` help topic and `start` on the `modelBayesControl()` help topic. No Bayesian calculations or dispatch behaviour change.

Stage 9.2.1 consumed build 1.1.1.004. Stage 9.2.2 therefore uses build 1.1.1.005.
