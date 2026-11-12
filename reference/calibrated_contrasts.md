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
# illness_death with a randomised treatment A0 added, for illustration.
set.seed(1)
treated <- illness_death
treated$A0 <- rbinom(max(treated$id), 1, 0.5)[treated$id]

fit <- fit_initial(
  treated,
  models = list(
    hazard_model(z ~ L0 + A0, limit = 1),
    death ~ L0 + A0 + z,
    censoring ~ L0
  ),
  tau = 1,
  treatment_model = A0 ~ L0
)

# \donttest{
# How much of the effect of treatment on death goes through its effect on
# the risk of z? The calibrated intervention gives the treated the risk of z
# they would have had under placebo (reference_a = 0):
contrasts <- calibrated_contrasts(
  treatment_fit = tmle_estimate(fit, risk(death), a = 1),
  placebo_fit = tmle_estimate(fit, risk(death), a = 0),
  calibrated_fit = tmle_estimate(
    fit,
    risk(death),
    calibrate(risk(z), relative_change(1), reference_a = 0),
    a = 1
  )
)
contrasts
#> Decomposition of the effect of a on death into mediated and calibrated contrasts
#> 
#>    target   contrast   estimate         se       lower      upper
#>    <char>     <char>      <num>      <num>       <num>      <num>
#> 1:  death      total -0.1405813 0.07724416 -0.29197704 0.01081449
#> 2:  death   mediated -0.0234265 0.01632300 -0.05541899 0.00856599
#> 3:  death calibrated -0.1171548 0.08141345 -0.27672220 0.04241266

# The total effect is the sum of the mediated and calibrated contrasts:
sum(contrasts$estimate$estimate[2:3]) - contrasts$estimate$estimate[1]
#> [1] 0
# }
```
