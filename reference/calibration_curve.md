# Trace the calibration curve over a grid of interventions

For each intervention `alpha` in `alpha_grid` on the intensity of the
process `z_name`, estimates the risk of `z_name` (the exposure) and of
each `target` under that intervention. Plotting the target estimates
against the exposure estimates traces out the calibration curve relating
the two. The standard errors of the target estimates account for
calibrating to the estimated exposure risk rather than to `alpha`.

## Usage

``` r
calibration_curve(
  initial_fit,
  target,
  z_name = "z",
  alpha_grid = seq(0, 5, length.out = 10),
  a = NULL,
  min_iter = 1,
  target_by_state = FALSE,
  block_size_z = NULL,
  block_size_target = NULL,
  cores = 1,
  verbose = FALSE,
  ...
)
```

## Arguments

- initial_fit:

  An initial fit, as returned by
  [`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md).

- target:

  Character vector. Names of the processes whose risks are estimated.

- z_name:

  String. Name of the process whose intensity `alpha` scales. If no such
  process exists in `initial_fit`, `alpha` is ignored with a warning.

- alpha_grid:

  Non-negative numeric vector. The interventions to evaluate.

- a:

  Optional baseline treatment value (`0` or `1`) to intervene on. Must
  be one of the values of `a` given to
  [`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md).

- min_iter:

  Non-negative integer. Minimum number of TMLE iterations, passed on to
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md).

- target_by_state:

  Logical. If `TRUE`, use state-specific clever covariates in the
  targeting step.

- block_size_z, block_size_target:

  Optional positive integers. If given, the expected time lost is
  estimated for `z_name` and `target`, respectively, instead of the
  risk; see `years_lost` in
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md).

- cores:

  Positive integer. Number of cores used by
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html) for
  the per-individual computations.

- verbose:

  Logical. If `TRUE`, print progress.

- ...:

  Further arguments passed on to
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md).

## Value

A list with elements:

- `estimate`: a `data.table` with one row per value of `alpha_grid` and
  columns `alpha`, `exposure_est`, `exposure_se`, `target_est` and
  `target_se` (suffixed by `_<target>` when there are several targets).

- `checks`: per value of `alpha_grid`, the estimated derivatives of the
  exposure and target estimates with respect to `alpha`.

- `eics`: per value of `alpha_grid`, the influence curves `alpha_eic`,
  `calibration_eic`, `target_eic` and `eic`.

## Examples

``` r
fit <- fit_initial(
  illness_death,
  fit_types = list(
    z = list(
      delta = 2,
      model = ~ L0,
      at_risk = function(dt) dt[["z"]] == 0
    ),
    death = list(delta = 1, model = ~ L0 + z),
    censoring = list(delta = 0, model = ~ L0)
  ),
  tau = 1
)

curve <- calibration_curve(
  fit,
  target = "death",
  z_name = "z",
  alpha_grid = c(0, 0.5, 1)
)
curve$estimate
#>    alpha exposure_est exposure_se target_est  target_se
#>    <num>        <num>       <num>      <num>      <num>
#> 1:   0.0   0.00000000  0.00000000  0.2820265 0.03975707
#> 2:   0.5   0.08791201  0.01791587  0.3040549 0.03896630
#> 3:   1.0   0.16448600  0.03131500  0.3239449 0.03949105
```
