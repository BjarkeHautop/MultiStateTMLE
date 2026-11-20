# Targeted estimation under a stochastic intensity intervention

Estimates a target parameter, such as the risk of death by the horizon
`tau` of the initial fit, under a stochastic intensity intervention,
optionally combined with an intervention `a` on the baseline treatment.
The initial fit from
[`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md)
is updated by iterative targeted maximum likelihood estimation (TMLE).

## Usage

``` r
tmle_estimate(
  initial_fit,
  target,
  intervention = NULL,
  a = NULL,
  control = list(),
  cores = 1,
  verbose = FALSE
)
```

## Arguments

- initial_fit:

  An initial fit, as returned by
  [`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md).

- target:

  The target parameter:
  [`risk()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/target-parameters.md)
  or
  [`years_lost()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/target-parameters.md)
  of one or more processes, each estimated separately under the same
  intervention.

- intervention:

  The intervention:
  [`scale_intensity()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/interventions.md),
  [`calibrate()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/interventions.md),
  or `NULL` (default) for none.

- a:

  Optional baseline treatment value (`0` or `1`) to intervene on. Must
  be one of the `arms` prepared by
  [`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md).
  Without it, each individual keeps their observed treatment.

- control:

  A named list of settings for how the estimate is computed:

  - `one_step` (`FALSE`): if `TRUE`, return the one-step estimator
    instead of iterating the TMLE update.

  - `block_size` (`10`): positive integer. For
    [`years_lost()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/target-parameters.md),
    the integral is approximated by evaluating the risk once per block
    of `block_size` consecutive time points. `1` gives the exact
    integral of the estimated risk; the computation time is roughly
    proportional to `1 / block_size`.

  - `target_by_state` (`FALSE`): if `TRUE`, use state-specific clever
    covariates in the targeting step.

  - `conv_const` (`1`): positive number. The TMLE iterations stop when
    the absolute mean of the efficient influence curve falls below
    `conv_const * se / log(n)`.

  - `max_iter` (`10`): maximum number of TMLE iterations.

  - `min_iter` (`0`): minimum number of TMLE iterations before checking
    convergence.

  - `truncate_weights` (`Inf`): positive number at which the clever
    weights of the intervention are truncated. `Inf` is no truncation.

- cores:

  Positive integer. Number of cores used by
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html) for
  the per-individual computations.

- verbose:

  Logical. If `TRUE`, print progress.

## Value

An object of class `tmle_estimate`, a list with elements:

- `estimate`: a `data.table` with one row per target process and columns
  `target`, `estimate`, its standard error `se`, and the bounds `lower`
  and `upper` of a 95% Wald confidence interval. Under a calibrated
  intervention, `se_crude` is the standard error that treats the scaling
  as known.

- `target`, `intervention` and `a`: as given.

- `method`: `"tmle"`, or `"one_step"` if `control$one_step = TRUE`.

- `eic`: a list, named by target process, with the efficient influence
  curve, one value per individual.

- `checks`: a list of diagnostics:

  - `iter`, `eic_solved_at` and `converged`: vectors, named by target
    process, of the number of iterations, the absolute mean of the
    efficient influence curve, and whether the iterations converged. A
    warning is given if they did not. Each iteration checks convergence
    before updating, so a converged fit made `iter - 1` targeting steps:
    `iter = 1` means the initial fit already met the criterion and was
    not updated.

  - `weights`: a matrix with the quantiles (columns `q50`, `q90`, `q99`
    and `q100`) over individuals of each individual's maximal
    `intervention`, `censoring` and `treatment` weight, and of their
    product `total`, the weight entering the targeting step (rows), as a
    positivity diagnostic. The intervention weights are those after
    truncation, and the treatment weights are 1 when `a` is not given.

  - `n_truncated`: the number of individuals whose initial intervention
    weight was truncated (0 if `control$truncate_weights = Inf`).

  `weights` and `n_truncated` are those of the first target process.

Without a calibrated intervention it also has `initial_estimate`, the
initial plug-in estimates before targeting, named by target process.
With one, it has:

- `calibration`: a `data.table` with the estimated scaling `alpha` and
  the exposure `theta` it attains, with standard errors and confidence
  intervals.

- `eic_alpha_known`: as `eic`, treating the scaling as known.

- `alpha_eic`: the efficient influence curve of the scaling.

- `checks$calibration`: a list with `converged` (whether a scaling
  attaining `theta` was found), `dist` (how far from `theta` it is),
  `c_n` (the tolerance), and the estimated derivatives with respect to
  the scaling of the exposure, `derivative_exposure`, and of the
  targets, `derivative_target`.

## Details

The intervention is either a fixed intensity scaling,
[`scale_intensity()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/interventions.md),
or a calibrated one,
[`calibrate()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/interventions.md),
whose scaling is chosen to bring an exposure to a given value. Under a
calibrated intervention the target is a composite parameter, and its
standard error accounts for the scaling being estimated.

## Examples

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

# The risk of death by tau = 1 without intervention, and if the intensity
# of illness were halved:
observed <- tmle_estimate(fit, risk(death))
halved <- tmle_estimate(fit, risk(death), scale_intensity(illness = 0.5))
observed
#> TMLE estimate of risk(death)
#> without intervention
#> 
#>    target  estimate         se     lower     upper
#>    <char>     <num>      <num>     <num>     <num>
#> 1:  death 0.3326759 0.03439736 0.2652584 0.4000935
#> 
#> Initial estimate: 0.3327; no targeting needed, as it already solves the EIC equation
halved
#> TMLE estimate of risk(death)
#> under scale_intensity(illness = 0.5)
#> 
#>    target  estimate         se     lower     upper
#>    <char>     <num>      <num>     <num>     <num>
#> 1:  death 0.2906665 0.03341164 0.2251809 0.3561521
#> 
#> Initial estimate: 0.2907; no targeting needed, as it already solves the EIC equation

# Several targets at once, here the risks of death and of illness themselves:
both <- tmle_estimate(
  fit,
  risk(death, illness),
  scale_intensity(illness = 0.5)
)
both$estimate
#>     target  estimate         se     lower     upper
#>     <char>     <num>      <num>     <num>     <num>
#> 1:   death 0.2906665 0.03341164 0.2251809 0.3561521
#> 2: illness 0.2169103 0.02428345 0.1693156 0.2645050

# The expected years lost before tau rather than the risk:
coef(tmle_estimate(fit, years_lost(death), scale_intensity(illness = 0.5)))
#>    death 
#> 0.149752 

# Diagnostics: the quantiles over individuals of the largest weights, as a
# check of positivity.
halved$checks$weights
#>                   q50      q90      q99     q100
#> intervention 1.142203 1.462249 1.906350 2.589352
#> censoring    1.110138 1.158804 1.201243 1.216979
#> treatment    1.000000 1.000000 1.000000 1.000000
#> total        1.249050 1.638933 2.228221 3.151186

# \donttest{
# A calibrated intervention: the risk of death if the risk of illness fell
# to 80% of its observed value. The first table gives the scaling alpha that
# achieves it and the risk of illness, theta, it attains:
cal <- tmle_estimate(
  fit,
  risk(death),
  calibrate(risk(illness), relative_change(0.8))
)
cal
#> TMLE estimate of risk(death)
#> under calibrate(risk(illness), relative_change(0.8))
#> 
#>    quantity  estimate         se     lower     upper
#>      <char>     <num>      <num>     <num>     <num>
#> 1:    alpha 0.7443751 0.01232345 0.7202216 0.7685287
#> 2:    theta 0.2950367 0.02790535 0.2403432 0.3497301
#> 
#> Estimates under the calibrated intervention:
#>    target  estimate         se    lower     upper   se_crude
#>    <char>     <num>      <num>    <num>     <num>      <num>
#> 1:  death 0.3127176 0.03330144 0.247448 0.3779872 0.03338458

# Its standard error accounts for alpha being estimated; se_crude treats
# alpha as known:
cal$estimate[, c("target", "estimate", "se", "se_crude")]
#>    target  estimate         se   se_crude
#>    <char>     <num>      <num>      <num>
#> 1:  death 0.3127176 0.03330144 0.03338458
# }
```
