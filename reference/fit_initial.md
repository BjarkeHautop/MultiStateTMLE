# Fit the initial hazard models

Fits the initial (nuisance) hazard model for each event process (Cox
models, or HAL), and the baseline treatment model if requested, and
expands the data to the grid of all observed event times up to `tau`. It
also builds the discrete product state space of the processes that other
processes' hazards depend on. The result is the input to
[`tmle_estimate()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_estimate.md)
and
[`intervention_curve()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/intervention_curve.md).

## Usage

``` r
fit_initial(
  data,
  models,
  tau,
  censoring = "censoring",
  treatment_model = NULL,
  arms = NULL,
  max_count = 1,
  control = list(),
  hal_control = list(),
  verbose = FALSE
)
```

## Arguments

- data:

  A `data.frame` of observed event history data in long format, with one
  row per event and columns `id`, `time` and `event` (the name of the
  process that jumped), plus any baseline covariates used in `models`.
  This is the format of
  [`simevent::simulate_events()`](https://rdrr.io/pkg/simevent/man/simulate_events.html).

- models:

  A list with the hazard model of each event process: a formula such as
  `death ~ biomarker + illness`, naming the process on the left and its
  covariates on the right, or a
  [`hazard_model()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/hazard_model.md)
  to also restrict when individuals are at risk of it. Other processes
  can be used as covariates by name, giving their number of events so
  far, and the time of an event through an indicator such as
  `I(event_time(illness, 1) > 0.5)` (see
  [`hazard_model()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/hazard_model.md)),
  also in the censoring model. Every event in `data` needs a model,
  including the censoring process.

- tau:

  Positive number. Follow-up horizon; the target parameters are
  evaluated at `tau`.

- censoring:

  String. Name of the censoring process in `models`.

- treatment_model:

  Optional logistic regression model for the baseline treatment, a
  two-sided formula such as `treatment ~ biomarker`. Required when
  `arms` is given.

- arms:

  Integer vector of the baseline treatment values (`0`/`1`) to prepare
  interventional hazards and clever weights for. Each value can then be
  passed as `a` to
  [`tmle_estimate()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_estimate.md)
  and
  [`intervention_curve()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/intervention_curve.md).
  Defaults to both values, `0:1`, when `treatment_model` is given.

- max_count:

  Positive integer. Maximum number of events tracked per recurrent
  process. Counts are capped at `max_count` in the state space, so
  hazards are assumed not to change after `max_count` events, and a
  recurrent process `z` in a model is fitted as the indicators
  `I(z >= 1) + ... + I(z >= max_count)`. It also bounds `k` in
  `event_time(z, k)`. A larger value is more flexible, but the state
  space, and with it the computation time, grows quickly with it.

- control:

  A named list of tuning settings:

  - `prune_states` (`FALSE`): if `TRUE`, prune states not admissible
    under the at-risk conditions of `models`.

  - `use_exponential` (`FALSE`): if `TRUE`, always transform the hazard
    increments `dH` to probabilities `1 - exp(-dH)`. Otherwise this is
    only done if some increments exceed 1.

- hal_control:

  A named list of settings for processes fitted with `fit = "hal"`, see
  Details.

- verbose:

  Logical. If `TRUE`, print model fits and progress.

## Value

An object of class `initial_fit`, a list with elements:

- `n`: the number of individuals.

- `tau`, `treatment_model` and `arms`: as given (`arms` with its default
  filled in).

- `models`: the hazard models, as a list of
  [`hazard_model()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/hazard_model.md)s
  named by process.

- `processes`: a `data.table` with one row per process and columns
  `name`, `type` (`"terminal"`, `"recurrent"` or `"one_jump"`, inferred
  from the data), `model` and `at_risk` (whether the process has a
  `limit` or `at_risk` condition).

- `cox_fits`: the fitted
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)
  model of each process, named by process. A recurrent process `z` in a
  model is fitted as the indicators
  `I(z >= 1) + ... + I(z >= max_count)`.
  [`as_sim_model()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/as_sim_model.md)
  builds a simulation model from them.

It also contains the internal elements used by
[`tmle_estimate()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_estimate.md),
among them `tmp_long` (the expanded data with fitted hazards and clever
weights) and `depend_matrix` (the product state space).

## Details

HAL fits are experimental and require a `fit.hal()` function to be
available on the search path. `hal_control` may contain `cut_time` (35),
`cut_one_way` (15), `cut_tk` (3), `cut_tk_values` (`NULL`), `two_way`
(`NULL`), `reduce_nk` (0.025), `sl` (a list of HAL library
specifications, built from `cut_one_way`, `cut_tk`, `cut_tk_values` and
`two_way` by default), `screen_two_way` (`FALSE`), `lambda_cvs`,
`event_dependent_cv` (`FALSE`), `npenalize_vars` (`NULL`), `v_folds`
(10), `seed` (`NULL`), `reduce_seed_dependence` (`FALSE`),
`penalize_time` (`FALSE`) and `verbose` (`FALSE`), with defaults in
parentheses.

Columns of `data` named after a process, such as the event counts that
[`simevent::simulate_events()`](https://rdrr.io/pkg/simevent/man/simulate_events.html)
adds, are ignored: the counts are recomputed from `event`.

## Parallel computation

The hazards in each state are predicted in parallel on the
[`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html) set
in the session, if any, such as `mirai::daemons(4)` for 4 local
processes.

## Examples

``` r
# Illness-death data: illness is an illness that can occur once, death ends
# follow-up, and so does censoring. One row per event:
head(illness_death)
#>   id      time     event  biomarker treatment
#> 1  1 2.8926838 censoring  0.2724785         1
#> 2  2 0.1282264   illness  0.3572619         0
#> 3  2 1.2558891     death  0.3572619         0
#> 4  3 0.6567113     death -0.8616620         0
#> 5  4 0.1752960   illness  0.8083350         1
#> 6  4 0.2409278     death  0.8083350         1

# One hazard model per process. The death hazard depends on whether the
# individual has fallen ill, and the illness can occur at most once:
fit <- fit_initial(
  illness_death,
  models = list(
    hazard_model(illness ~ biomarker + treatment, limit = 1),
    death ~ biomarker + treatment + illness,
    censoring ~ biomarker
  ),
  tau = 1
)
fit
#> Initial fit of 200 individuals up to tau = 1
#> 
#>         name     type                                   model at_risk
#>       <char>   <char>                                  <char>  <lgcl>
#> 1:   illness one_jump         illness ~ biomarker + treatment    TRUE
#> 2:     death terminal death ~ biomarker + treatment + illness   FALSE
#> 3: censoring terminal                   censoring ~ biomarker   FALSE

# The process types are inferred from the data, and the hazard models are
# kept for reference:
fit$processes
#>         name     type                                   model at_risk
#>       <char>   <char>                                  <char>  <lgcl>
#> 1:   illness one_jump         illness ~ biomarker + treatment    TRUE
#> 2:     death terminal death ~ biomarker + treatment + illness   FALSE
#> 3: censoring terminal                   censoring ~ biomarker   FALSE
fit$models$death
#> Hazard model: death ~ biomarker + treatment + illness

# With a baseline treatment, a treatment model prepares estimates under
# both arms, for `a` in tmle_estimate():
fit_treated <- fit_initial(
  illness_death,
  models = list(
    hazard_model(illness ~ biomarker + treatment, limit = 1),
    death ~ biomarker + treatment + illness,
    censoring ~ biomarker
  ),
  tau = 1,
  treatment_model = treatment ~ biomarker
)
fit_treated
#> Initial fit of 200 individuals up to tau = 1
#> 
#>         name     type                                   model at_risk
#>       <char>   <char>                                  <char>  <lgcl>
#> 1:   illness one_jump         illness ~ biomarker + treatment    TRUE
#> 2:     death terminal death ~ biomarker + treatment + illness   FALSE
#> 3: censoring terminal                   censoring ~ biomarker   FALSE
#> 
#> Treatment model: treatment ~ biomarker, prepared for arms a = 0, 1

# \donttest{
# Data simulated with simevent can be passed as is. Here relapse is
# recurrent, and the death hazard grows with the number of relapses so far,
# tracked up to max_count:
if (requireNamespace("simevent", quietly = TRUE)) {
  model <- simevent::sim_model(
    biomarker = simevent::sim_covariate(function(N) rnorm(N)),
    relapse = simevent::sim_process("transient", eta = 0.6, nu = 1),
    death = simevent::sim_process("terminal", eta = 0.2, nu = 1),
    censoring = simevent::sim_process("censoring", eta = 0.2, nu = 1),
    effects = list(
      simevent::sim_effect(relapse ~ biomarker, coef = 0.5),
      simevent::sim_effect(death ~ relapse, coef = 0.4)
    )
  )
  relapses <- simevent::simulate_events(model, n = 200, seed = 1)

  fit_initial(
    relapses,
    models = list(
      relapse ~ biomarker,
      death ~ biomarker + relapse,
      censoring ~ biomarker
    ),
    tau = 2,
    max_count = 3
  )
}
#> Initial fit of 200 individuals up to tau = 2
#> 
#>         name      type                       model at_risk
#>       <char>    <char>                      <char>  <lgcl>
#> 1:   relapse recurrent         relapse ~ biomarker   FALSE
#> 2:     death  terminal death ~ biomarker + relapse   FALSE
#> 3: censoring  terminal       censoring ~ biomarker   FALSE
# }
```
