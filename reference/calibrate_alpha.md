# Calibrate an intensity intervention to a target value

Finds the intervention `alpha` on the intensity of the process `z_name`
such that the estimated risk of `z_name` equals a calibration target,
and estimates the risk of `target` under that calibrated intervention.
The standard error of the latter accounts for the estimation of `alpha`.

## Usage

``` r
calibrate_alpha(
  initial_fit,
  target,
  z_name = "z",
  theta = NULL,
  rho = NULL,
  delta = NULL,
  rho_1a = NULL,
  a = NULL,
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

  Character vector. Names of the processes whose risks are estimated
  under the calibrated intervention.

- z_name:

  String. Name of the process whose intensity `alpha` scales. If no such
  process exists in `initial_fit`, `alpha` is ignored with a warning.

- theta, rho, delta, rho_1a:

  Numbers specifying the calibration target, see Description. Give
  exactly one.

- a:

  Optional baseline treatment value (`0` or `1`) to intervene on. Must
  be one of the values of `a` given to
  [`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md).

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

- `estimate`: a named vector with the calibrated `alpha_est` and its
  `alpha_se`, the risk of `target` (`target_est`) with standard error
  `target_se` (accounting for the estimation of alpha) and
  `target_se_crude` (treating alpha as known), and the calibration
  target `theta` with standard error `theta_se`.

- `checks`: diagnostics of the alpha search and derivative estimates.

- `alpha_eic`, `target_eic` and `eic`: the influence curves of the
  calibrated alpha, of the risk of `target` treating alpha as known, and
  of the risk of `target`.

With several targets, the per-target entries of `estimate` and `checks`
are suffixed by the target name (e.g. `target_est_death`), and
`target_eic` and `eic` are lists named by target.

## Details

Exactly one of `theta`, `rho`, `delta` and `rho_1a` specifies the
calibration target:

- `theta`: the risk of `z_name` itself.

- `rho`: a fraction of the risk of `z_name` without intervention.

- `delta`: an absolute difference from the risk of `z_name` without
  intervention.

- `rho_1a`: a fraction of the risk of `z_name` without intervention on
  `z_name`, under the opposite baseline treatment `1 - a`.

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

# Risk of death under the intervention on z that brings the risk of z
# down to 80% of its observed value.
cal <- calibrate_alpha(fit, target = "death", z_name = "z", rho = 0.8)
cal$estimate
#>       alpha_est        alpha_se      target_est       target_se target_se_crude 
#>      0.79272380      0.00995552      0.31992034      0.03864198      0.03866085 
#>           theta        theta_se 
#>      0.12713026      0.02505712 
```
