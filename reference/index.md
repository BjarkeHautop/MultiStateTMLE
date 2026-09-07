# Package index

## Package overview

Overview of the MultiStateTMLE package.

- [`MultiStateTMLE`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/MultiStateTMLE-package.md)
  [`MultiStateTMLE-package`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/MultiStateTMLE-package.md)
  : MultiStateTMLE: Causal Inference for Multistate Event Processes via
  TMLE

## Simulating data

Functions for simulating multistate event history data, either from a
user-specified process/effect specification, or from parameters
previously fitted to observed data.

- [`sim.generic()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/sim.generic.md)
  : Simulate event history data from a generic process specification
- [`sim.from.data()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/sim.from.data.md)
  : Simulate event history data from parameters fitted to observed data

## Initial fit

Fitting the initial (nuisance) hazard models and building the discrete
product state space used internally.

- [`prepare.initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/prepare.initial.md)
  : Prepare the initial fit for multistate TMLE
- [`compute.Q.clever.per.id()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/compute.Q.clever.per.id.md)
  : Compute Q and clever covariates per id using a discrete-state
  backward recursion

## Targeted estimation

Targeting the initial fit under a stochastic intensity intervention.

- [`tmle.alpha.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle.alpha.fun.md)
  : Targeted maximum likelihood estimation under a shape-parameter
  intervention
- [`estimate.derivative()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/estimate.derivative.md)
  : Estimate the numerical derivative of an alpha-indexed estimand

## Calibration

Calibrating a shape-parameter intervention to a target effect size, and
decomposing treatment effects into calibrated contrasts.

- [`estimate.alpha.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/estimate.alpha.fun.md)
  : Search for the alpha achieving a target estimand value
- [`calibration.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibration.fun.md)
  : Calibrate a shape-parameter intervention to a target effect size
- [`calibration.curve.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibration.curve.fun.md)
  : Compute a calibration curve over a grid of shape-parameter
  interventions
- [`make.calibrated.contrasts()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/make.calibrated.contrasts.md)
  : Compute total, mediated, and calibrated treatment contrasts
