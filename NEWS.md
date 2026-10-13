# MultiStateTMLE (development version)

## Deprecations

- All exported functions have snake_case replacements, with snake_case
  arguments. The old functions now warn via `lifecycle::deprecate_warn()` and
  will be removed in a future release:
  - `prepare.initial()` → `fit_initial()`
  - `tmle.alpha.fun()` → `tmle_alpha()`
  - `calibration.fun()` → `calibrate_alpha()`
  - `calibration.curve.fun()` → `calibration_curve()`
  - `estimate.alpha.fun()` → `solve_alpha()`
  - `estimate.derivative()` → `estimate_derivative()`
  - `make.calibrated.contrasts()` → `calibrated_contrasts()`
  - `compute.Q.clever.per.id()` is now internal.

## Breaking changes

- Removed `sim.generic()` and `sim.from.data()`, which duplicated simevent's
  functions of the same name. Use `simevent::sim_model()`,
  `simevent::sim_events()` and `simevent::sim_model_from_fits()` instead.
  Accordingly, `prepare.initial(return.parameters.for.simulation = TRUE)` is no
  longer supported.

## Changes in the new functions

Compared with the functions they replace, the new functions:

- Validate their arguments, and default to `cores = 1` and `verbose = FALSE`.
  Progress is reported with `message()` and only when `verbose = TRUE`.
- No longer have `browse` arguments.
- `fit_initial()` takes `data` and `fit_types` first and requires `tau` and
  `fit_types`. It renames `fit.treatment` to `treatment_model` and accepts
  `at_risk` in `fit_types`. It collects the HAL settings into a `hal_control`
  list, drops `use.cores.prediction`, `max.Tk` and `cv.glmnet`, and folds
  `verbose.exponential` into `verbose`. Model strings may contain whitespace,
  and the `delta` code is no longer limited to one digit.
- `tmle_alpha()` requires `initial_fit` and `target`. It no longer fits the
  initial model itself, so the `tau` and `dt` arguments are gone; the horizon is
  always that of `initial_fit` (as it already was when an initial fit was
  supplied). The same applies to `tau` and `tau.z` in `calibrate_alpha()` and
  `calibration_curve()`, which had no effect.
- `calibrate_alpha()` requires exactly one of `theta`, `rho`, `delta` and
  `rho_1a` instead of defaulting to `theta = 0.5`.
- `estimate_derivative()` takes `fun` first, has no defaults for `fun` and `h`,
  and returns an unnamed number. Consequently the `checks` of
  `calibrate_alpha()` are named `est.deriv.auxiliary`/`est.deriv.target`, and
  `theta` is no longer suffixed with `.tmle.est`.
- `solve_alpha()` drops the unused `use.cores` argument and now respects
  `trace_every`.

## Bug fixes

- Fixed `mclapply()` erroring with a negative number of cores on machines with 5
  or fewer cores.
- The influence curve at `alpha = 0` in `calibration_curve()` (and at
  `theta = 0` in `calibrate_alpha()`) is no longer `NaN`.

## New features

- Added the simulated example dataset `illness_death`.

# MultiStateTMLE 0.0.0.9000

- Initial development version.
