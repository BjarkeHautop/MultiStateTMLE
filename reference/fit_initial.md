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
  arms = NULL,
  max_count = 1,
  derived_vars = list(),
  depend_time = list(),
  control = list(),
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

  A named list with one entry per event process, named by the process.
  Each entry is a list with elements:

  - `delta`: the process's integer code in the `delta` column of `data`.
    Exactly one process must have code `0` (censoring).

  - `model`: a one-sided formula of the covariates of its hazard model,
    such as `~ L0 + z`. Other processes can be used as covariates by
    name, in which case their counting process (lagged by one event) is
    used.

  - `fit`: `"cox"` (default) or `"hal"`.

  - `at_risk`: an optional function of the data returning a logical
    vector, marking the rows at risk of the process (e.g.
    `function(dt) dt[["z"]] == 0` for a process that can jump at most
    once).

- tau:

  Positive number. Follow-up horizon; the target parameters are
  evaluated at `tau`.

- treatment_model:

  Optional logistic regression model for the baseline treatment, a
  two-sided formula such as `A0 ~ L0`. Required when `arms` is given.

- arms:

  Integer vector of the baseline treatment values (`0`/`1`) to prepare
  interventional hazards and clever weights for. Each value can then be
  passed as `a` to
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md),
  [`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md)
  and
  [`calibration_curve()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibration_curve.md).
  Defaults to both values, `0:1`, when `treatment_model` is given.

- max_count:

  Positive integer. Maximum number of jumps tracked per recurrent
  process in the state space.

- derived_vars:

  Named list of functions of the data, each deriving an additional state
  variable, named by the list, that can be used in the `model`s of
  `fit_types`.

- depend_time:

  Experimental, and may change. Named list of functions making hazards
  depend on the time of a jump. Each is named by a state variable: a
  process `z` that jumps once, or `z.k` for the `k`th jump of a
  recurrent process `z`. It adds the covariate `T.<name>`, `0` before
  that jump and the function applied to the time of the jump after it,
  which can be used (as a factor) in the `model`s of `fit_types`.

- control:

  A named list of tuning settings:

  - `prune_states` (`FALSE`): if `TRUE`, prune states not admissible
    under the `at_risk` functions of `fit_types`.

  - `use_exponential` (`FALSE`): if `TRUE`, always transform the hazard
    increments `dH` to probabilities `1 - exp(-dH)`. Otherwise this is
    only done if some increments exceed 1.

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

An object of class `initial_fit`, a list with elements:

- `n`: the number of individuals.

- `tau`, `treatment_model` and `arms`: as given (`arms` with its default
  filled in).

- `processes`: a `data.table` with one row per process and columns
  `name`, `delta`, `type` (`"terminal"`, `"recurrent"` or `"one_jump"`,
  inferred from the data), `model` and `at_risk` (whether an `at_risk`
  function was given).

It also contains the internal elements used by
[`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md),
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

## Functions in the arguments

The functions given to `fit_initial()` and
[`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md)
are called in different ways:

- `at_risk` in `fit_types`, and the functions in `derived_vars`, are
  called with a `data.table` and must return one value per row. They are
  called both with the data, with a row per individual and time, and
  with the state space, with a row per state and a column per process
  that other processes depend on. So they may only use those process
  columns, such as `dt[["z"]]`, the number of jumps of `z` so far.

- The functions in `depend_time` are called with a numeric vector of
  times and must return one value per time.

- `target_only_in_state` in
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md)
  is called with the state space only, `initial_fit$depend_matrix`, and
  must return one logical per state.

- The functions in `alpha_list` in
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md)
  are called with the columns of the expanded data
  `initial_fit$tmp_long` named by their arguments, such as `time`, a
  baseline covariate or a process, and must return one value per row.

## Examples

``` r
# Illness-death data: z is an illness that can occur once, death is
# terminal, and delta = 0 is censoring.
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
      delta = 2,
      model = ~ L0,
      at_risk = function(dt) dt[["z"]] == 0
    ),
    death = list(delta = 1, model = ~ L0 + z),
    censoring = list(delta = 0, model = ~ L0)
  ),
  tau = 1
)
fit
#> Initial fit of 150 individuals up to tau = 1
#> 
#>         name delta     type   model at_risk
#>       <char> <num>   <char>  <char>  <lgcl>
#> 1:         z     2 one_jump     ~L0    TRUE
#> 2:     death     1 terminal ~L0 + z   FALSE
#> 3: censoring     0 terminal     ~L0   FALSE
```
