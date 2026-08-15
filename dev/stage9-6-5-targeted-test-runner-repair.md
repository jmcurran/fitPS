# Stage 9.6.5 targeted test runner repair

Stage 9.6.4 correctly stabilized the external Poisson-normal likelihood, but its targeted validation step called `testthat::test_file()` directly. That does not load the fitPS package or the normal test helper environment, so the targeted test failed immediately because package functions such as `makePSData()` and helper-defined model constructors were unavailable.

Stage 9.6.5 changes only the stage validation workflow. The targeted regression now uses `devtools::test(filter = "generic-numerical-posterior", ...)`, which loads fitPS and its test helpers in the same environment as the full package test suite while restricting execution to the generic numerical posterior regression file. The full test suite still runs only after the targeted regression passes.

No Bayesian mathematics, cubature implementation, model contract, or user-facing package behavior is changed by this repair.
