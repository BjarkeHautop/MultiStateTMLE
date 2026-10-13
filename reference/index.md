# Package index

## Package overview

Overview of the MultiStateTMLE package.

- [`MultiStateTMLE`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/MultiStateTMLE-package.md)
  [`MultiStateTMLE-package`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/MultiStateTMLE-package.md)
  : MultiStateTMLE: Causal Inference for Multistate Event Processes via
  TMLE

## Initial fit

Fitting the initial (nuisance) hazard models and building the discrete
product state space used internally.

- [`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md)
  : Fit the initial hazard models

## Targeted estimation

Targeting the initial fit under a stochastic intensity intervention.

- [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md)
  : Targeted estimation under a stochastic intensity intervention
- [`estimate_derivative()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/estimate_derivative.md)
  : Estimate the derivative of an estimate with respect to alpha

## Calibration

Calibrating an intensity intervention to a target effect size, and
decomposing treatment effects into calibrated contrasts.

- [`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md)
  : Calibrate an intensity intervention to a target value
- [`calibration_curve()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibration_curve.md)
  : Trace the calibration curve over a grid of interventions
- [`solve_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/solve_alpha.md)
  : Solve for the intervention achieving a target value
- [`calibrated_contrasts()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrated_contrasts.md)
  : Decompose a treatment effect into mediated and calibrated contrasts

## Data

Example data.

- [`illness_death`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/illness_death.md)
  : Simulated illness-death data

## Deprecated

Functions superseded by the snake_case functions above. Will be removed
in the future.

- [`prepare.initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/prepare.initial.md)
  : Prepare the initial fit (deprecated)
- [`tmle.alpha.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle.alpha.fun.md)
  : TMLE under a shape-parameter intervention (deprecated)
- [`calibration.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibration.fun.md)
  : Calibrate a shape-parameter intervention (deprecated)
- [`calibration.curve.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibration.curve.fun.md)
  : Calibration curve over a grid of interventions (deprecated)
- [`estimate.alpha.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/estimate.alpha.fun.md)
  : Search for the alpha achieving a target value (deprecated)
- [`estimate.derivative()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/estimate.derivative.md)
  : Estimate the derivative with respect to alpha (deprecated)
- [`make.calibrated.contrasts()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/make.calibrated.contrasts.md)
  : Total, mediated and calibrated contrasts (deprecated)
- [`compute.Q.clever.per.id()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/compute.Q.clever.per.id.md)
  : Compute Q and clever covariates per id (deprecated)
