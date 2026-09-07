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
