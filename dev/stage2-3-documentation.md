# Stage 2.3 documentation audit

Stage 2.3 follows the Stage 2.2 parameterisation refactor by aligning user-facing documentation with the package convention:

- `shape` is the standard zeta parameter alpha.
- Valid zeta shape values satisfy `shape > 1`.
- VGAM uses a shifted internal parameter equal to `shape - 1`.
- fitPS converts to the VGAM shifted parameter only at VGAM-backed zeta function boundaries.

The main documentation repairs are in `fitDist()`, `fitZIDist()`, `makePrior()`, README, and NEWS. The obsolete statement that printed shape differed from stored shape should not reappear; fitted objects now store and print standard shape values.
