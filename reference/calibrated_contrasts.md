# Decompose a treatment effect into mediated and calibrated contrasts

Given estimates of the risk of an outcome under treatment, under
placebo, and under a calibrated intervention (as returned by
[`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md)
and
[`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md)),
computes the total treatment effect and its decomposition into a
mediated contrast (treatment versus calibrated) and a calibrated
contrast (calibrated versus placebo), with Wald confidence intervals
based on the influence curves.

## Usage

``` r
calibrated_contrasts(
  treatment_fit,
  placebo_fit,
  calibrated_fit,
  target = NULL,
  conf_level = 0.95
)
```

## Arguments

- treatment_fit, placebo_fit:

  Fits for the treatment and placebo arms, as returned by
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md).

- calibrated_fit:

  Fit for the calibrated intervention, as returned by
  [`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md).

- target:

  Optional string. If `calibrated_fit` contains estimates for several
  targets, the name of the one to use.

- conf_level:

  Number between 0 and 1. Confidence level of the intervals.

## Value

A `data.table` with one row per contrast (`total`, `mediated` and
`calibrated`) and columns `target` (`NA` if `target` is `NULL`),
`contrast`, `estimate`, `se`, `lower` and `upper`. The attributes
`decomposition_error` and `eic_decomposition_error` give how far the
estimates and influence curves are from satisfying total = mediated +
calibrated, which should be zero up to numerical error.

## Examples

``` r
set.seed(1405)
treatment_fit <- list(estimate = c(tmle_est = 0.7), eic = rnorm(200, sd = 0.1))
placebo_fit <- list(estimate = c(tmle_est = 0.3), eic = rnorm(200, sd = 0.1))
calibrated_fit <- list(estimate = c(target_est = 0.5), eic = rnorm(200, sd = 0.1))

calibrated_contrasts(treatment_fit, placebo_fit, calibrated_fit)
#>    target   contrast estimate         se     lower     upper
#>    <char>     <char>    <num>      <num>     <num>     <num>
#> 1:   <NA>      total      0.4 0.01028048 0.3798506 0.4201494
#> 2:   <NA>   mediated      0.2 0.01050568 0.1794092 0.2205908
#> 3:   <NA> calibrated      0.2 0.01028145 0.1798487 0.2201513
```
