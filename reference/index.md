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
- [`hazard_model()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/hazard_model.md)
  : Specify the hazard model of an event process
- [`as_sim_model()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/as_sim_model.md)
  : Build a simulation model from an initial fit

## Targeted estimation

Estimating a target parameter under a stochastic intensity intervention,
fixed or calibrated.

- [`tmle_estimate()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_estimate.md)
  : Targeted estimation under a stochastic intensity intervention
- [`risk()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/target-parameters.md)
  [`years_lost()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/target-parameters.md)
  : Target parameters
- [`scale_intensity()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/interventions.md)
  [`calibrate()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/interventions.md)
  [`level()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/interventions.md)
  [`absolute_change()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/interventions.md)
  [`relative_change()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/interventions.md)
  : Stochastic intensity interventions

## Calibration

Tracing the intervention curve, and decomposing treatment effects into
calibrated contrasts.

- [`intervention_curve()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/intervention_curve.md)
  : Trace the intervention curve over a grid of intensity scalings
- [`calibrated_contrasts()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrated_contrasts.md)
  : Decompose a treatment effect into mediated and calibrated contrasts

## Results

Printing, extracting estimates from, and computing confidence intervals
for the result objects.

- [`print(`*`<initial_fit>`*`)`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/result-methods.md)
  [`print(`*`<tmle_estimate>`*`)`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/result-methods.md)
  [`print(`*`<intervention_curve>`*`)`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/result-methods.md)
  [`print(`*`<calibrated_contrasts>`*`)`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/result-methods.md)
  [`coef(`*`<tmle_estimate>`*`)`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/result-methods.md)
  [`confint(`*`<tmle_estimate>`*`)`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/result-methods.md)
  : Methods for the result objects

## Data

Example data.

- [`illness_death`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/illness_death.md)
  : Simulated illness-death trial data
