# Stochastic intensity interventions

Describe the intervention that
[`tmle_estimate()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_estimate.md)
estimates the target under.

## Usage

``` r
scale_intensity(..., only_first = FALSE)

calibrate(exposure, toward, reference_a = NULL)

level(theta)

absolute_change(delta)

relative_change(rho)
```

## Arguments

- ...:

  Named factors, one per process to intervene on, such as `z = 0.5`.
  Each is a non-negative number, or a one-sided formula giving a
  time-varying factor, such as `z = ~ ifelse(time < 0.5, 0.5, 1)`, that
  may use `time`, baseline covariates and the processes (their number of
  events so far).

- only_first:

  Logical. If `TRUE`, follow-up ends at the second event of an
  intervened process, so that only its first event is intervened on.

- exposure:

  The exposure to calibrate, a
  [`risk()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/target-parameters.md)
  or
  [`years_lost()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/target-parameters.md)
  of the single process whose intensity is scaled.

- toward:

  The calibration rule: `level()`, `absolute_change()` or
  `relative_change()`.

- reference_a:

  Optional baseline treatment value (`0` or `1`) under which the
  exposure without intervention is estimated, for `absolute_change()`
  and `relative_change()`. Defaults to the `a` given to
  [`tmle_estimate()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_estimate.md).
  A different value calibrates across arms, e.g. the risk of the
  exposure under treatment to that under placebo.

- theta:

  Number between 0 and 1. The exposure to attain.

- delta:

  Number. The change in the exposure.

- rho:

  Non-negative number. The factor on the exposure.

## Value

`scale_intensity()` and `calibrate()` return an object of class
`tmle_intervention`; `level()`, `absolute_change()` and
`relative_change()` one of class `calibration_rule`.

## Details

`scale_intensity()` is a fixed intensity-scaling intervention: it
multiplies the intensity of each named process by its factor `alpha`,
\\\Lambda^z \mapsto \alpha \Lambda^z\\. `alpha = 1` is no intervention,
`alpha < 1` attenuates the process and `alpha > 1` intensifies it.

`calibrate()` is a calibrated intervention: the intensity of the process
of `exposure` is scaled by the `alpha` that brings the exposure to the
value given by `toward`. The target is then a composite parameter, whose
standard error accounts for `alpha` being estimated. The calibration
rules are:

- `level(theta)`: the exposure equals `theta`.

- `absolute_change(delta)`: the exposure is `delta` plus its value
  without intervention.

- `relative_change(rho)`: the exposure is `rho` times its value without
  intervention.

## Examples

``` r
# Fixed intensity scalings. Halve the intensity of z, remove z altogether,
# or double it:
scale_intensity(z = 0.5)
#> Intervention: scale_intensity(z = 0.5)
scale_intensity(z = 0)
#> Intervention: scale_intensity(z = 0)
scale_intensity(z = 2)
#> Intervention: scale_intensity(z = 2)

# A time-varying scaling: halve the intensity of z before time 0.5 only.
scale_intensity(z = ~ ifelse(time < 0.5, 0.5, 1))
#> Intervention: scale_intensity(z = ~ifelse(time < 0.5, 0.5, 1))

# A scaling that depends on a baseline covariate: halve it only for
# individuals with L0 > 0.
scale_intensity(z = ~ ifelse(L0 > 0, 0.5, 1))
#> Intervention: scale_intensity(z = ~ifelse(L0 > 0, 0.5, 1))

# Several processes at once:
scale_intensity(z = 0.5, death = 0.8)
#> Intervention: scale_intensity(z = 0.5, death = 0.8)

# Calibrated interventions, scaling the intensity of z so that its risk
# equals 0.1, falls by 0.05, or falls to 80% of its value without
# intervention:
calibrate(risk(z), level(0.1))
#> Intervention: calibrate(risk(z), level(0.1))
calibrate(risk(z), absolute_change(-0.05))
#> Intervention: calibrate(risk(z), absolute_change(-0.05))
calibrate(risk(z), relative_change(0.8))
#> Intervention: calibrate(risk(z), relative_change(0.8))

# Use them in tmle_estimate():
fit <- fit_initial(
  illness_death,
  models = list(
    hazard_model(z ~ L0, limit = 1),
    death ~ L0 + z,
    censoring ~ L0
  ),
  tau = 1
)
tmle_estimate(fit, risk(death), scale_intensity(z = ~ ifelse(L0 > 0, 0.5, 1)))
#> TMLE estimate of risk(death)
#> under scale_intensity(z = ~ifelse(L0 > 0, 0.5, 1))
#> 
#>    target  estimate        se     lower     upper
#>    <char>     <num>     <num>     <num>     <num>
#> 1:  death 0.3159039 0.0390286 0.2394092 0.3923985
#> 
#> Initial estimate: 0.3159; converged after 1 iterations
# \donttest{
tmle_estimate(fit, risk(death), calibrate(risk(z), level(0.1)))
#> TMLE estimate of risk(death)
#> under calibrate(risk(z), level(0.1))
#> 
#>    quantity  estimate       se     lower     upper
#>      <char>     <num>    <num>     <num>     <num>
#> 1:    alpha 0.6157885 0.139007 0.3433398 0.8882371
#> 2:    theta 0.1000000 0.000000 0.1000000 0.1000000
#> 
#> Estimates under the calibrated intervention:
#>    target estimate         se     lower     upper   se_crude
#>    <char>    <num>      <num>     <num>     <num>      <num>
#> 1:  death 0.315415 0.03871762 0.2395298 0.3913001 0.03850806
# }
```
