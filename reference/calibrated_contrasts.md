# Decompose a treatment effect into mediated and calibrated contrasts

Given estimates of the risk of an outcome under treatment, under
placebo, and under a calibrated intervention (as returned by
[`tmle_estimate()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_estimate.md)),
computes the total treatment effect and its decomposition into a
mediated contrast (treatment versus calibrated) and a calibrated
contrast (calibrated versus placebo), with Wald confidence intervals
based on the influence curves.

## Usage

``` r
calibrated_contrasts(treatment_fit, placebo_fit, calibrated_fit, level = 0.95)
```

## Arguments

- treatment_fit, placebo_fit:

  Fits for the treatment and placebo arms, as returned by
  [`tmle_estimate()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_estimate.md)
  with different baseline treatments `a` and no calibrated intervention.

- calibrated_fit:

  Fit under a
  [`calibrate()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/interventions.md)
  intervention, as returned by
  [`tmle_estimate()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_estimate.md)
  with the `a` of `treatment_fit`. It may contain estimates for several
  target processes; the one of `treatment_fit` is used.

- level:

  Number between 0 and 1. Confidence level of the intervals.

## Value

An object of class `calibrated_contrasts`, a list with elements:

- `estimate`: a `data.table` with one row per contrast (`total`,
  `mediated` and `calibrated`) and columns `target`, `contrast`,
  `estimate`, `se`, `lower` and `upper`.

- `target` and `level`: the target and the confidence level.

- `checks`: a list with `decomposition_error` and
  `eic_decomposition_error`, how far the estimates and influence curves
  are from satisfying total = mediated + calibrated, which should be
  zero up to numerical error.

## Examples

``` r
# A trial where treatment lowers the risk of illness, which raises
# the risk of death:
fit <- fit_initial(
  illness_death,
  models = list(
    hazard_model(illness ~ biomarker + treatment, limit = 1),
    death ~ biomarker + treatment + illness,
    censoring ~ biomarker
  ),
  tau = 1,
  treatment_model = treatment ~ biomarker
)

# \donttest{
# How much of the effect of treatment on death goes through its effect on
# the risk of illness? The calibrated intervention gives the treated the
# risk of illness they would have had under placebo (reference_a = 0):
contrasts <- calibrated_contrasts(
  treatment_fit = tmle_estimate(fit, risk(death), a = 1),
  placebo_fit = tmle_estimate(fit, risk(death), a = 0),
  calibrated_fit = tmle_estimate(
    fit,
    risk(death),
    calibrate(risk(illness), relative_change(1), reference_a = 0),
    a = 1
  )
)
contrasts
#> Decomposition of the effect of a on death into mediated and calibrated contrasts
#> 
#>    target   contrast    estimate         se      lower         upper
#>    <char>     <char>       <num>      <num>      <num>         <num>
#> 1:  death      total -0.18273318 0.06913193 -0.3182293 -0.0472370913
#> 2:  death   mediated -0.06387838 0.03275458 -0.1280762  0.0003194037
#> 3:  death calibrated -0.11885480 0.08063644 -0.2768993  0.0391897229

# The total effect is the sum of the mediated and calibrated contrasts:
sum(contrasts$estimate$estimate[2:3]) - contrasts$estimate$estimate[1]
#> [1] 0
# }
```
