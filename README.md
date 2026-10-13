
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
- `tmle_alpha()` targets the initial fit under a stochastic intensity
  intervention, which scales the intensity of a process by `alpha`
  (optionally combined with an intervention on baseline treatment). It
  returns a TMLE estimate, its standard error, and optionally the
  efficient influence curve.
- `calibrate_alpha()` and `calibration_curve()` calibrate the
  intervention to a target effect size, and `calibrated_contrasts()`
  decomposes a total treatment effect into mediated and calibrated
  contrasts.

For example, with the simulated illness-death data `illness_death`,
where an illness `z` raises the risk of death `outcome1`, the risk of
death by time 1 if the intensity of the illness were halved is estimated
by:

``` r
fit <- fit_initial(
  illness_death,
  fit_types = list(
    z = list(
      model = "Surv(tstart, tstop, delta == 2) ~ L0",
      at_risk = function(dt) dt[["z"]] == 0
    ),
    outcome1 = "Surv(tstart, tstop, delta == 1) ~ L0 + z",
    censoring = "Surv(tstart, tstop, delta == 0) ~ L0"
  ),
  tau = 1
)

est <- tmle_alpha(fit, target = "outcome1", alpha = 0.5, z_name = "z")
est$estimate
#>   tmle.est         se      g.est 
#> 0.30405495 0.03853505 0.31239982
```

To simulate multistate event history data, see the
[simevent](https://github.com/BjarkeHautop/simevent) package.

## Contributing

All contributions are welcome!
