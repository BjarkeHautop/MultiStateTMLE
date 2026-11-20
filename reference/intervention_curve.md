# Trace the intervention curve over a grid of intensity scalings

For each scaling `alpha` of the intensity of the exposure process,
estimates the exposure and the target under
[`scale_intensity()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/interventions.md)
with that `alpha`. Plotting the target estimates against the exposure
estimates traces out how the target responds to changes in the exposure
on its own scale, such as its risk. The standard errors of the target
estimates account for calibrating to the estimated exposure rather than
to `alpha`, as in
[`calibrate()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/interventions.md).

## Usage

``` r
intervention_curve(
  initial_fit,
  target,
  exposure,
  alpha = seq(0, 5, length.out = 10),
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
  of one or more processes.

- exposure:

  The exposure:
  [`risk()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/target-parameters.md)
  or
  [`years_lost()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/target-parameters.md)
  of the single process whose intensity is scaled.

- alpha:

  Non-negative numeric vector. The scalings to evaluate.

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

An object of class `intervention_curve`, a list with elements:

- `estimate`: a `data.table` with one row per value of `alpha` and
  target process, and columns `alpha`, `target`, the exposure
  (`exposure_estimate`, with standard error `exposure_se`), and the
  target (`estimate`, with standard error `se` and the bounds `lower`
  and `upper` of a 95% Wald confidence interval).

- `target`, `exposure` and `a`: as given.

- `checks`: a `data.table` with the estimated derivatives with respect
  to `alpha` of the exposure (`derivative_exposure`) and target
  (`derivative_target`) estimates, per `alpha` and target process.

- `eic`: a list with, per value of `alpha`, the influence curves of the
  exposure estimate (`exposure_eic`), of `alpha` calibrated to it
  (`alpha_eic`), and lists named by target process of the target
  estimates accounting for its calibration (`eic`) and treating `alpha`
  as known (`eic_alpha_known`), as in
  [`tmle_estimate()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_estimate.md).

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

# \donttest{
# How the risk of death responds to the risk of illness, as the intensity
# of illness is scaled from 0 (no illness) to 2 (twice the observed
# intensity):
curve <- intervention_curve(
  fit,
  risk(death),
  exposure = risk(illness),
  alpha = c(0, 0.5, 1, 2)
)
curve
#> Intervention curve of risk(death) against risk(illness)
#> over 4 scalings of the intensity of illness
#> 
#>    alpha target exposure_estimate exposure_se  estimate         se     lower
#>    <num> <char>             <num>       <num>     <num>      <num>     <num>
#> 1:   0.0  death         0.0000000  0.00000000 0.2336750 0.03966975 0.1559237
#> 2:   0.5  death         0.2169103  0.02428345 0.2906665 0.03374676 0.2245241
#> 3:   1.0  death         0.3687958  0.03488168 0.3326759 0.03427321 0.2655017
#> 4:   2.0  death         0.5616300  0.04251611 0.3893948 0.04055059 0.3099171
#>        upper
#>        <num>
#> 1: 0.3114263
#> 2: 0.3568089
#> 3: 0.3998502
#> 4: 0.4688724

# Plot the risk of death against the risk of illness, with 95% confidence
# intervals:
est <- curve$estimate
plot(
  est$exposure_estimate,
  est$estimate,
  type = "b",
  ylim = range(est$lower, est$upper),
  xlab = "Risk of illness by time 1",
  ylab = "Risk of death by time 1"
)
arrows(
  est$exposure_estimate,
  est$lower,
  est$exposure_estimate,
  est$upper,
  angle = 90,
  code = 3,
  length = 0.03
)

# }
```
