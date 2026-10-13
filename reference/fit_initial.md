# Fit the initial hazard models

Fits the initial (nuisance) hazard model for each event process (Cox
models, or HAL), and the baseline treatment model if requested, and
expands the data to the grid of all observed event times up to `tau`. It
also builds the discrete product state space of the processes that other
processes' hazards depend on. The result is the input to
[`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md),
[`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md)
and
[`calibration_curve()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibration_curve.md).

## Usage

``` r
fit_initial(
  data,
  fit_types,
  tau,
  treatment_model = NULL,
  a = NULL,
  max_count = 1,
  derived_vars = list(),
  depend_time = list(),
  prune_states = FALSE,
  use_exponential = FALSE,
  hal_control = list(),
  cores = 1,
  verbose = FALSE
)
```

## Arguments

- data:

  A `data.frame` of observed event history data in long format, with one
  row per event and columns `id`, `time` and `delta` (the integer code
  of the process that jumped, with `0` reserved for censoring), plus any
  baseline covariates used in `fit_types`.

- fit_types:

  A named list with one entry per event process, each either a model
  formula string or a list with elements:

  - `model`: a formula string of the form
    `"Surv(tstart, tstop, delta == k) ~ covariates"`, where `k` is the
    process's `delta` code. Other processes can be used as covariates by
    name, in which case their counting process (lagged by one event) is
    used.

  - `fit`: `"cox"` (default) or `"hal"`.

  - `at_risk`: an optional function of the data returning a logical
    vector, marking the rows at risk of the process (e.g.
    `function(dt) dt[["z"]] == 0` for a process that can jump at most
    once).

  Exactly one process must have `delta` code `0` (censoring).

- tau:

  Positive number. Follow-up horizon; the target parameters are
  evaluated at `tau`.

- treatment_model:

  Optional model for the baseline treatment `A0`, either a formula
  string such as `"A0 ~ L0"` or a list with elements `model` and `fit`
  (only `"glm"` is supported). Required when `a` is given.

- a:

  Optional integer vector of baseline treatment values (`0`/`1`) to
  compute interventional hazards and clever weights for. Each value can
  then be passed as `a` to
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md).

- max_count:

  Positive integer. Maximum number of jumps tracked per recurrent
  process in the state space.

- derived_vars:

  Named list of functions of the data deriving additional state
  variables.

- depend_time:

  Named list of functions describing dependence on the time since a jump
  of a process, keyed by process name.

- prune_states:

  Logical. If `TRUE`, prune states not admissible under the `at_risk`
  functions of `fit_types`.

- use_exponential:

  Logical. If `TRUE`, always transform the hazard increments `dH` to
  probabilities `1 - exp(-dH)`. Otherwise this is only done if some
  increments exceed 1.

- hal_control:

  A named list of settings for processes fitted with `fit = "hal"`, see
  Details.

- cores:

  Positive integer. Number of cores used by
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  when predicting hazards across states.

- verbose:

  Logical. If `TRUE`, print model fits and progress.

## Value

A list containing, among others, `tmp.long` (the expanded data with
fitted hazards and clever weights), `depend.matrix` (the product state
space), `process.names`, `process.types`, `process.deltas`,
`cens.process.id` and `at.risks`.

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

## Examples

``` r
# Illness-death data: z is an illness that can occur once, outcome1 is
# death, and delta = 0 is censoring.
head(illness_death)
#>   id      time delta         L0
#> 1  1 1.0836890     1  0.2724785
#> 2  2 0.6539194     1  0.3572619
#> 3  3 0.1212156     1 -0.8616620
#> 4  4 1.1416233     2  0.8083350
#> 5  4 3.6778032     1  0.8083350
#> 6  5 1.0962284     2  0.6127352

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
fit$process.types
#> $z
#> [1] "one.jump"
#> 
#> $outcome1
#> [1] "terminal"
#> 
#> $censoring
#> [1] "terminal"
#> 
```
