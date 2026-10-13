# CRAN compliance

- HAL fits call `fit.hal()`, which is not part of the package.

# Cleanup

- Internal variable names are still dotted.

- [`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md)
  parses process `delta` codes and covariates out of model strings;
  consider specifying processes by name, matching simevent’s
  `sim_events()` output (`event` column).

# Performance:

- R/fit-initial.R: repeated
  [`grep()`](https://rdrr.io/r/base/grep.html)/[`setdiff()`](https://rdrr.io/r/base/sets.html)/[`unique()`](https://rdrr.io/r/base/unique.html)
  recomputation on the same vectors inside inner loops. Needs a closer
  look with tests before touching.
