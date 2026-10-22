# Targeted estimation under a stochastic intensity intervention

Estimates a target parameter under a stochastic intensity intervention
that scales the intensity of the process `z_name` by `alpha`, optionally
combined with an intervention `a` on the baseline treatment. The initial
fit from
[`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md)
is updated by iterative targeted maximum likelihood estimation (TMLE).

## Usage

``` r
tmle_alpha(
  initial_fit,
  target,
  z_name = NULL,
  alpha = 1,
  alpha_list = NULL,
  a = NULL,
  estimand = c("risk", "years_lost"),
  only_first = FALSE,
  target_only_in_state = NULL,
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

  Character vector. Names of the target processes.

- z_name:

  String. Name of the process in `initial_fit` whose intensity `alpha`
  scales. Required unless `alpha = 1` or `alpha_list` is given.

- alpha:

  Non-negative number. Multiplicative intervention on the intensity of
  `z_name`. `alpha = 1` is no intervention.

- alpha_list:

  Experimental, and may change. Optional named list of functions, one
  per process, giving a time-varying multiplicative intervention on that
  process's intensity. Each function is called with the columns of
  `initial_fit$tmp_long` named by its arguments, e.g. `time`, a baseline
  covariate such as `L0`, or a process such as `z` (its number of jumps
  so far), and must return one non-negative factor per row. For example,
  `list(z = function(time) ifelse(time < 0.5, 0.5, 1))` halves the
  intensity of `z` before time 0.5. Cannot be combined with `alpha` and
  `z_name`.

- a:

  Optional baseline treatment value (`0` or `1`) to intervene on. Must
  be one of the `arms` prepared by
  [`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md).

- estimand:

  The target parameter: `"risk"`, the risk of `target` by `tau` (or the
  expected number of events for a recurrent `target`), or
  `"years_lost"`, the expected time lost, i.e. the integral of the risk
  of `target` up to `tau`.

- only_first:

  Logical. If `TRUE`, only the first `target` event counts, and only the
  first jump of `z_name` is intervened on.

- target_only_in_state:

  Experimental, and may change. Optional function of the state space
  `initial_fit$depend_matrix`, a data frame with one row per state and a
  column per process that other processes depend on (its number of
  jumps), returning a logical vector marking the states in which
  `target` events count. For example,
  `function(states) states[["z"]] == 0` counts only deaths before
  illness.

- control:

  A named list of settings for how the estimate is computed:

  - `one_step` (`FALSE`): if `TRUE`, return the one-step estimator
    instead of iterating the TMLE update.

  - `block_size` (`10`): positive integer. For
    `estimand = "years_lost"`, the integral is approximated by
    evaluating the risk once per block of `block_size` consecutive time
    points. `1` gives the exact integral of the estimated risk; the
    computation time is roughly proportional to `1 / block_size`.

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

An object of class `tmle_alpha`, a list with elements:

- `estimate`: a `data.table` with one row per `target` and columns
  `target`, `estimate`, its standard error `se`, and the bounds `lower`
  and `upper` of a 95% Wald confidence interval.

- `initial_estimate`: the initial plug-in estimates, before targeting,
  named by target.

- `method`: `"tmle"`, or `"one_step"` if `control$one_step = TRUE`.

- `estimand`: as given.

- `alpha`, `z_name` and `a`: the intervention (`alpha` is `NULL` if
  `alpha_list` is given).

- `eic`: a list, named by `target`, with the efficient influence curve,
  one value per individual.

- `checks`: a list of diagnostics:

  - `iter`, `eic_solved_at` and `converged`: vectors, named by target,
    of the number of iterations, the absolute mean of the efficient
    influence curve, and whether the iterations converged. A warning is
    given if they did not.

  - `weights`: a matrix with the quantiles (columns `q50`, `q90`, `q99`
    and `q100`) over individuals of each individual's maximal
    `intervention`, `censoring` and `treatment` weight, and of their
    product `total`, the weight entering the targeting step (rows), as a
    positivity diagnostic. The intervention weights are those after
    truncation, and the treatment weights are 1 when `a` is not given.

  - `n_truncated`: the number of individuals whose initial intervention
    weight was truncated (0 if `control$truncate_weights = Inf`).

  `weights` and `n_truncated` are shared by the targets; with
  `only_first` and several targets, they are those of the first.

## Details

By default the target parameter is the risk of a `target` event by the
horizon `tau` given to
[`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md)
(or the expected number of events for a recurrent `target`). Several
targets are estimated separately, under the same intervention.

The functions in `alpha_list` and `target_only_in_state` are called as
described in the section "Functions in the arguments" of
[`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md).

## Examples

``` r
fit <- fit_initial(
  illness_death,
  fit_types = list(
    z = list(
      delta = 2,
      model = ~ L0,
      at_risk = function(dt) dt[["z"]] == 0
    ),
    death = list(delta = 1, model = ~ L0 + z),
    censoring = list(delta = 0, model = ~ L0)
  ),
  tau = 1
)

# Risk of death by tau = 1 if the intensity of z were halved.
est <- tmle_alpha(fit, target = "death", z_name = "z", alpha = 0.5)
est$estimate
#>    target  estimate         se     lower     upper
#>    <char>     <num>      <num>     <num>     <num>
#> 1:  death 0.3040549 0.03853505 0.2285276 0.3795823
```
