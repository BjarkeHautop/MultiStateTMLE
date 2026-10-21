# MultiStateTMLE

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

- [`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md)
  fits the initial (nuisance) hazard models for each event process in
  the data, and builds the discrete product state space used internally.
- [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md)
  targets the initial fit under a stochastic intensity intervention,
  which scales the intensity of a process by `alpha` (optionally
  combined with an intervention on baseline treatment). It returns a
  TMLE estimate, its standard error, and optionally the efficient
  influence curve.
- [`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md)
  and
  [`calibration_curve()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibration_curve.md)
  calibrate the intervention to a target effect size, and
  [`calibrated_contrasts()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrated_contrasts.md)
  decomposes a total treatment effect into mediated and calibrated
  contrasts.

For example, with the simulated illness-death data `illness_death`,
where an illness `z` raises the risk of `death`, the risk of death by
time 1 if the intensity of the illness were halved is estimated by:

``` r

fit <- fit_initial(
  illness_death,
  fit_types = list(
    z = list(
      delta = 2,
      model = ~L0,
      at_risk = function(dt) dt[["z"]] == 0
    ),
    death = list(delta = 1, model = ~ L0 + z),
    censoring = list(delta = 0, model = ~L0)
  ),
  tau = 1
)

est <- tmle_alpha(fit, target = "death", alpha = 0.5, z_name = "z")
est$estimate
#>   tmle.est         se      g.est 
#> 0.30405495 0.03853505 0.31239982
```

To simulate multistate event history data, see the
[simevent](https://github.com/BjarkeHautop/simevent) package.

## Contributing

All contributions are welcome!
