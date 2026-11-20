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

  Named factors, one per process to intervene on, such as
  `illness = 0.5`. Each is a non-negative number, or a one-sided formula
  giving a time-varying factor, such as
  `illness = ~ ifelse(time < 0.5, 0.5, 1)`, that may use `time`,
  baseline covariates and the processes (their number of events so far).

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
# Fixed intensity scalings. Halve the intensity of illness, remove illness
# altogether, or double it:
scale_intensity(illness = 0.5)
#> Intervention: scale_intensity(illness = 0.5)
scale_intensity(illness = 0)
#> Intervention: scale_intensity(illness = 0)
scale_intensity(illness = 2)
#> Intervention: scale_intensity(illness = 2)

# A time-varying scaling: halve the intensity of illness before time 0.5
# only.
scale_intensity(illness = ~ ifelse(time < 0.5, 0.5, 1))
#> Intervention: scale_intensity(illness = ~ifelse(time < 0.5, 0.5, 1))

# A scaling that depends on a baseline covariate: halve it only for
# individuals with biomarker > 0.
scale_intensity(illness = ~ ifelse(biomarker > 0, 0.5, 1))
#> Intervention: scale_intensity(illness = ~ifelse(biomarker > 0, 0.5, 1))

# Several processes at once:
scale_intensity(illness = 0.5, death = 0.8)
#> Intervention: scale_intensity(illness = 0.5, death = 0.8)

# Calibrated interventions, scaling the intensity of illness so that its risk
# equals 0.1, falls by 0.05, or falls to 80% of its value without
# intervention:
calibrate(risk(illness), level(0.1))
#> Intervention: calibrate(risk(illness), level(0.1))
calibrate(risk(illness), absolute_change(-0.05))
#> Intervention: calibrate(risk(illness), absolute_change(-0.05))
calibrate(risk(illness), relative_change(0.8))
#> Intervention: calibrate(risk(illness), relative_change(0.8))

# Use them in tmle_estimate():
fit <- fit_initial(
  illness_death,
  models = list(
    hazard_model(illness ~ biomarker + treatment, limit = 1),
    death ~ biomarker + treatment + illness,
    censoring ~ biomarker
  ),
  tau = 1
)
tmle_estimate(
  fit,
  risk(death),
  scale_intensity(illness = ~ ifelse(biomarker > 0, 0.5, 1))
)
#> TMLE estimate of risk(death)
#> under scale_intensity(illness = ~ifelse(biomarker > 0, 0.5, 1))
#> 
#>    target  estimate         se     lower     upper
#>    <char>     <num>      <num>     <num>     <num>
#> 1:  death 0.3053054 0.03359683 0.2394568 0.3711539
#> 
#> Initial estimate: 0.3053; no targeting needed, as it already solves the EIC equation
# \donttest{
tmle_estimate(fit, risk(death), calibrate(risk(illness), level(0.1)))
#> TMLE estimate of risk(death)
#> under calibrate(risk(illness), level(0.1))
#> 
#>    quantity  estimate         se     lower     upper
#>      <char>     <num>      <num>     <num>     <num>
#> 1:    alpha 0.2176164 0.02987141 0.1590696 0.2761633
#> 2:    theta 0.1000000 0.00000000 0.1000000 0.1000000
#> 
#> Estimates under the calibrated intervention:
#>    target  estimate         se    lower     upper   se_crude
#>    <char>     <num>      <num>    <num>     <num>      <num>
#> 1:  death 0.2607342 0.03592421 0.190324 0.3311443 0.03559196
# }
```
