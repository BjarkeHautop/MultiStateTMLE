# Compute total, mediated, and calibrated treatment contrasts

Given TMLE fits for a treatment arm, a placebo arm, and a calibrated
(intermediate) intervention, computes the total treatment effect and its
decomposition into a mediated and a calibrated contrast, with influence
function-based confidence intervals.

## Usage

``` r
make.calibrated.contrasts(
  target = NULL,
  treatment.fit,
  placebo.fit,
  calibrated.fit,
  conf.level = 0.95
)
```

## Arguments

- target:

  optional character; name of a specific target within `calibrated.fit`
  to use (when `calibrated.fit` contains estimates/EICs for multiple
  targets).

- treatment.fit:

  TMLE fit object for the treatment arm, as returned by
  [`tmle.alpha.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle.alpha.fun.md).

- placebo.fit:

  TMLE fit object for the placebo/control arm.

- calibrated.fit:

  TMLE fit object for the calibrated intervention.

- conf.level:

  confidence level for the Wald-type confidence intervals.

## Value

A `data.table` with one row per contrast (`total`, `mediated`,
`calibrated`) and columns `estimate`, `se`, `lower`, `upper`, with
`decomposition.error` and `eic.decomposition.error` attributes.

## Examples

``` r
set.seed(1405)
treatment.fit <- list(estimate = c(tmle.est = 0.7), eic = rnorm(200, sd = 0.1))
placebo.fit <- list(estimate = c(tmle.est = 0.3), eic = rnorm(200, sd = 0.1))
calibrated.fit <- list(estimate = c(target.est = 0.5), eic = rnorm(200, sd = 0.1))

make.calibrated.contrasts(
  treatment.fit = treatment.fit,
  placebo.fit = placebo.fit,
  calibrated.fit = calibrated.fit
)
#>      contrast estimate         se     lower     upper
#>        <char>    <num>      <num>     <num>     <num>
#> 1:      total      0.4 0.01028048 0.3798506 0.4201494
#> 2:   mediated      0.2 0.01050568 0.1794092 0.2205908
#> 3: calibrated      0.2 0.01028145 0.1798487 0.2201513
```
