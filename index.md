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

- [`prepare.initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/prepare.initial.md)
  fits the initial (nuisance) hazard models for each process in the data
  (Cox, HAL, or exponential), and builds the discrete product state
  space used internally.
- [`tmle.alpha.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle.alpha.fun.md)
  targets the initial fit under a stochastic intensity intervention
  (scaling the shape of a process by `alpha`, optionally combined with
  an intervention on baseline treatment), returning a TMLE estimate, its
  standard error, and (optionally) the efficient influence curve.
- [`calibration.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibration.fun.md)/[`calibration.curve.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibration.curve.fun.md)
  calibrate the shape intervention to a target effect size, and
  [`make.calibrated.contrasts()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/make.calibrated.contrasts.md)
  decomposes a total treatment effect into mediated and calibrated
  contrasts.

[`sim.generic()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/sim.generic.md)
and
[`sim.from.data()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/sim.from.data.md)
(built on top of [simevent](https://github.com/BjarkeHautop/simevent))
simulate multistate event history data, either from a user-specified
process/effect specification, or from parameters previously fitted to
observed data.

## Contributing

All contributions are welcome!
