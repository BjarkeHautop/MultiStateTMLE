
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MultiStateTMLE

<!-- badges: start -->

[![R-CMD-check](https://github.com/BjarkeHautop/MultiStateTMLE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BjarkeHautop/MultiStateTMLE/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/BjarkeHautop/MultiStateTMLE/graph/badge.svg)](https://app.codecov.io/gh/BjarkeHautop/MultiStateTMLE)
<!-- badges: end -->

MultiStateTMLE is a framework for nonparametric causal inference with
multistate event processes, built around stochastic intensity
interventions, calibrated interventions, and composite estimands, and
estimated using targeted maximum likelihood estimation (TMLE) allowing
flexible nuisance learning.

## Installation

You can install the development version of MultiStateTMLE from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("BjarkeHautop/MultiStateTMLE")
```

## Usage

``` r
library(MultiStateTMLE)
```

The pipeline has three main steps:

- `fit_initial()` fits the initial (nuisance) hazard models for each
  event process in the data, and builds the discrete product state space
  used internally.
- `tmle_estimate()` targets the initial fit to estimate a target, such
  as `risk(death)`, under a stochastic intensity intervention: a fixed
  scaling of the intensity of a process,
  `scale_intensity(illness = 0.5)`, or one calibrated to bring an
  exposure to a given value,
  `calibrate(risk(illness), relative_change(0.8))` (optionally combined
  with an intervention on baseline treatment). It returns a TMLE
  estimate with its standard error and efficient influence curve.
- `intervention_curve()` traces how the target responds to the exposure
  over a grid of scalings, and `calibrated_contrasts()` decomposes a
  total treatment effect into mediated and calibrated contrasts.

For example, with the simulated illness-death data `illness_death`,
where `illness` raises the risk of `death`, the risk of death by time 1
if the intensity of the illness were halved is estimated by:

``` r
fit <- fit_initial(
  illness_death,
  models = list(
    hazard_model(illness ~ biomarker + treatment, limit = 1),
    death ~ biomarker + treatment + illness,
    censoring ~ biomarker
  ),
  tau = 1
)

tmle_estimate(fit, risk(death), scale_intensity(illness = 0.5))
#> TMLE estimate of risk(death)
#> under scale_intensity(illness = 0.5)
#>
#>    target  estimate         se     lower     upper
#>    <char>     <num>      <num>     <num>     <num>
#> 1:  death 0.2906665 0.03341164 0.2251809 0.3561521
#>
#> Initial estimate: 0.2907; no targeting needed, as it already solves the EIC equation
```

To simulate multistate event history data, see the
[simevent](https://github.com/BjarkeHautop/simevent) package, whose
output `fit_initial()` takes as is.

## Contributing

All contributions are welcome!
