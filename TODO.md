# CRAN compliance

- HAL fits call `fit.hal()`, which is not part of the package.

# Cleanup

- Internal variable names are still dotted.

- `fit_initial()` parses process `delta` codes and covariates out of model
  strings; consider specifying processes by name, matching simevent's
  `sim_events()` output (`event` column).

# Performance:

- R/fit-initial.R: repeated `grep()`/`setdiff()`/`unique()` recomputation on the
  same vectors inside inner loops. Needs a closer look with tests before
  touching.
