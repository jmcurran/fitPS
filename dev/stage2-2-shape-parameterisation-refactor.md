# Stage 2.2 shape parameterisation refactor notes

Stage 2.2 changes fitPS so package-facing `shape` means the standard zeta parameter alpha, with `shape > 1`.

VGAM still uses a shifted parameter internally. fitPS now converts at VGAM boundaries with:

```r
vgamShape = shape - 1
```

The active tests in `tests/testthat/test-zeta-shape-parameterisation.R` check that probability functions, prediction, fitted values, and random generation use the standard shape convention.
