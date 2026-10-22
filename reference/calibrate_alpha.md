# Calibrate an intensity intervention to a target value

Finds the intervention `alpha` on the intensity of the process `z_name`
such that the estimated risk of `z_name` (the exposure) equals a
calibration target `theta`, and estimates the risk of each `target`
under that calibrated intervention. The standard errors of the latter
account for the estimation of `alpha`.

## Usage

``` r
calibrate_alpha(
  initial_fit,
  target,
  z_name,
  value,
  scale = c("absolute", "ratio", "difference"),
  a = NULL,
  reference_a = a,
  estimand = c("risk", "years_lost"),
  exposure_estimand = c("risk", "years_lost"),
  control = list(),
  cores = 1,
  verbose = FALSE
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

  String. Name of the process in `initial_fit` whose intensity is
  intervened on and calibrated.

- value:

  Number. The calibration target on the `scale`, see Description. Must
  be non-negative unless `scale = "difference"`.

- scale:

  The scale of `value`: `"absolute"`, `"ratio"` or `"difference"`, see
  Description.

- a:

  Optional baseline treatment value (`0` or `1`) to intervene on, as in
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md).

- reference_a:

  Optional baseline treatment value (`0` or `1`) under which the
  reference exposure of a `"ratio"` or `"difference"` `scale` is
  estimated. Defaults to `a`; e.g. `a = 1, reference_a = 0` calibrates
  the exposure under treatment to that under placebo.

- estimand:

  The target parameter for `target`, as in
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md).

- exposure_estimand:

  The target parameter for the exposure `z_name`, that `value` and the
  calibration refer to, as in
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md).

- control:

  Settings for the estimates, passed on to
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md).

- cores:

  Positive integer. Number of cores used by
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html) for
  the per-individual computations.

- verbose:

  Logical. If `TRUE`, print progress.

## Value

An object of class `calibrate_alpha`, a list with elements:

- `estimate`: a `data.table` with one row per `target` and columns
  `target`, the risk `estimate` under the calibrated intervention, its
  standard error `se` (accounting for the estimation of alpha), the
  bounds `lower` and `upper` of a 95% Wald confidence interval, and the
  standard error `se_crude` treating alpha as known.

- `calibration`: a `data.table` with the calibrated `alpha` and the
  calibration target `theta` (column `quantity`), with `estimate`, `se`,
  `lower` and `upper`.

- `z_name`, `a`, `value`, `scale`, `reference_a`, `estimand` and
  `exposure_estimand`: as given (`reference_a` is `NULL` for the
  `"absolute"` `scale`).

- `eic`, `eic_alpha_known`: lists, named by target, of the influence
  curves of the risk of `target`, accounting for the estimation of alpha
  and treating it as known.

- `alpha_eic`: the influence curve of the calibrated alpha.

- `checks`: diagnostics of the alpha search (`converged`, `dist`, `c_n`)
  and the estimated derivatives with respect to alpha of the risk of
  `z_name` (`derivative_z`) and of each target (`derivative_target`).

## Details

The calibration target `theta` is given by `value` on the `scale`:

- `"absolute"`: `theta = value`.

- `"ratio"`: `theta = value * reference`.

- `"difference"`: `theta = value + reference`.

where `reference` is the exposure without intervention on `z_name`,
under the baseline treatment `reference_a`.

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
cal <- calibrate_alpha(
  fit,
  target = "death",
  z_name = "z",
  value = 0.8,
  scale = "ratio"
)
cal$estimate
#>    target  estimate         se     lower     upper   se_crude
#>    <char>     <num>      <num>     <num>     <num>      <num>
#> 1:  death 0.3199203 0.03864198 0.2441835 0.3956572 0.03866085
```
