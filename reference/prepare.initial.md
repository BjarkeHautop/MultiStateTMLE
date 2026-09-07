# Prepare the initial fit for multistate TMLE

Fits initial (nuisance) hazard models for each process in `fit.types`
(Cox, HAL, or exponential), builds the discrete product state space and
its dependency structure, and computes the baseline clever
weights/predictions needed by
[`tmle.alpha.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle.alpha.fun.md).
This is the initialization step of the multistate TMLE pipeline.

## Usage

``` r
prepare.initial(
  dt,
  tau = 1.2,
  fit.types = list(Z = list(model = "Surv(tstart, tstop, delta == 2)~L0+L+A0", fit =
    "cox", at.risk = function(dt) (dt[["Z"]] == 0)), L = list(model =
    "Surv(tstart, tstop, delta == 3)~L0+Z+A0", fit = "cox", at.risk = function(dt)
    (dt[["L"]] == 0)), outcome1 = list(model =
    "Surv(tstart, tstop, delta == 1)~L0+Z+L+A0", fit = "cox"), censoring = list(model =
    "Surv(tstart, tstop, delta == 0)~L0+Z+L+A0", fit = "cox")),
  use.exponential = FALSE,
  verbose.exponential = FALSE,
  browse = FALSE,
  max.count = 1,
  a = NULL,
  derived.vars = list(),
  depend.time = list(),
  fit.treatment = NULL,
  prune.states = FALSE,
  cut.time = 35,
  cut.one.way = 15,
  cut.Tk = 3,
  cut.Tk.values = NULL,
  max.Tk = max.count,
  two.way = NULL,
  reduce.NK = 0.025,
  hal.sl = list(c(cut.one.way = cut.one.way, cut.Tk = cut.Tk, cut.Tk.values =
    cut.Tk.values, two.way = two.way)),
  screen.two.way = FALSE,
  lambda.cvs = c((9:1)/10, (9:2)/10^2, seq(1/10^2, 1/10^5, length = 100)),
  event.dependent.cv = FALSE,
  npenalize.vars = NULL,
  V = 10,
  seed.hal = NULL,
  reduce.seed.dependence = FALSE,
  penalize.time = FALSE,
  use.cores = 50,
  use.cores.prediction = 1,
  verbose.hal = FALSE,
  browse.hal = FALSE,
  cv.glmnet = FALSE,
  verbose = FALSE,
  return.parameters.for.simulation = FALSE
)
```

## Arguments

- dt:

  data.table of observed event history data (long format, with `id`,
  `time`, `delta`, and one column per process).

- tau:

  follow-up horizon (in time units).

- fit.types:

  named list describing, for each process, the hazard `model` (a
  `Surv(...)` formula string), the `fit` method (`"cox"`, `"hal"`, or
  `"glm"`), and an optional `at.risk` function.

- use.exponential:

  logical; if `TRUE`, use exponential (rather than Cox/HAL) hazard fits
  where applicable.

- verbose.exponential:

  logical; if `TRUE`, print diagnostics for the exponential fits.

- browse:

  logical; if `TRUE`, drop into
  [`browser()`](https://rdrr.io/r/base/browser.html).

- max.count:

  maximum recurrent-event count tracked per process in the product state
  space.

- a:

  optional intervention on the baseline treatment `A0`; if supplied,
  predictions are also computed under `A0` set to each value in `a`.

- derived.vars:

  list of functions deriving auxiliary state variables from the
  discovered processes.

- depend.time:

  list describing additional time-dependency structure between
  processes, used when HAL is not fit.

- fit.treatment:

  optional model specification (as for `fit.types`) for the baseline
  treatment `A0` propensity model.

- prune.states:

  logical; if `TRUE`, prune unreachable/zero-probability states from the
  product state space.

- cut.time:

  HAL basis cut count for the time variable.

- cut.one.way:

  HAL basis cut count for one-way (main effect) terms.

- cut.Tk:

  HAL basis cut count for history/time-since-event terms.

- cut.Tk.values:

  optional explicit cut values for `cut.Tk`.

- max.Tk:

  maximum history order for HAL history terms.

- two.way:

  optional list of two-way HAL interaction term specifications.

- reduce.NK:

  numeric; quantile used to cap the recurrent-event count `NK` when
  reducing HAL basis dimension.

- hal.sl:

  list of HAL specifications, one per Super Learner library candidate.

- screen.two.way:

  logical; if `TRUE`, screen two-way HAL terms before fitting.

- lambda.cvs:

  numeric vector of lasso penalty values to cross-validate over when
  fitting HAL.

- event.dependent.cv:

  logical; if `TRUE`, stratify HAL cross-validation folds by event
  status.

- npenalize.vars:

  optional character vector of HAL basis terms to exempt from lasso
  penalization.

- V:

  number of cross-validation folds for HAL.

- seed.hal:

  optional random seed for HAL cross-validation.

- reduce.seed.dependence:

  logical; if `TRUE`, average over multiple HAL cross-validation seeds.

- penalize.time:

  logical; if `TRUE`, include the time variable itself among the
  penalized HAL basis terms.

- use.cores:

  number of cores to use for model fitting.

- use.cores.prediction:

  number of cores to use for prediction steps.

- verbose.hal:

  logical; if `TRUE`, print HAL fitting diagnostics.

- browse.hal:

  logical; if `TRUE`, drop into
  [`browser()`](https://rdrr.io/r/base/browser.html) inside the HAL
  fitting step.

- cv.glmnet:

  logical; currently not supported.

- verbose:

  logical; if `TRUE`, print overall progress.

- return.parameters.for.simulation:

  logical; if `TRUE`, additionally return fitted parameters formatted
  for
  [`sim.from.data()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/sim.from.data.md).

## Value

A list (the `initial.fit` object) containing, among others, `tmp.long`
(the long-format data with fitted hazards/weights), `depend.matrix` (the
product state space), `process.names`, `process.types`,
`process.deltas`, `cens.process.id`, and `at.risks`.
