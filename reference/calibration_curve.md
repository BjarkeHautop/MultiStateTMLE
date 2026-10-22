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
  z_name,
  alpha_grid = seq(0, 5, length.out = 10),
  a = NULL,
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

  Character vector. Names of the processes whose risks are estimated.

- z_name:

  String. Name of the process in `initial_fit` whose intensity
  `alpha_grid` scales.

- alpha_grid:

  Non-negative numeric vector. The interventions to evaluate.

- a:

  Optional baseline treatment value (`0` or `1`) to intervene on. Must
  be one of the `arms` prepared by
  [`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md).

- estimand:

  The target parameter for `target`, as in
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md).

- exposure_estimand:

  The target parameter for the exposure `z_name`, as in
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

An object of class `calibration_curve`, a list with elements:

- `estimate`: a `data.table` with one row per value of `alpha_grid` and
  `target`, and columns `alpha`, `target`, the risk of `z_name`
  (`exposure_estimate`, with standard error `exposure_se`), and the risk
  of `target` (`estimate`, with standard error `se` and the bounds
  `lower` and `upper` of a 95% Wald confidence interval).

- `z_name`, `a`, `estimand` and `exposure_estimand`: as given.

- `checks`: a `data.table` with the estimated derivatives with respect
  to `alpha` of the exposure (`derivative_exposure`) and target
  (`derivative_target`) estimates, per `alpha` and `target`.

- `eic`: a list with, per value of `alpha_grid`, the influence curves of
  the exposure estimate (`exposure_eic`), of `alpha` calibrated to it
  (`alpha_eic`), and lists named by target of the target estimates
  accounting for its calibration (`eic`) and treating `alpha` as known
  (`eic_alpha_known`), as in
  [`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md).

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
#>    alpha target exposure_estimate exposure_se  estimate         se     lower
#>    <num> <char>             <num>       <num>     <num>      <num>     <num>
#> 1:   0.0  death        0.00000000  0.00000000 0.2820265 0.03975707 0.2041041
#> 2:   0.5  death        0.08486244  0.01787691 0.3040549 0.03899060 0.2276348
#> 3:   1.0  death        0.15891283  0.03132140 0.3250469 0.03915660 0.2483014
#>        upper
#>        <num>
#> 1: 0.3599489
#> 2: 0.3804751
#> 3: 0.4017924
```
