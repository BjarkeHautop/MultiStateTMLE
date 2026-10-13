# Prepare the initial fit (deprecated)

`prepare.initial()` is deprecated as of MultiStateTMLE 0.1.0. Use
[`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md)
instead.

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
  use.cores = 1,
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

  Observed data, see `data` in
  [`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md).

- tau, a, verbose:

  See
  [`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md).

- fit.types, max.count, derived.vars, depend.time, prune.states,
  use.exponential:

  See the snake_case equivalents in
  [`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md).
  Use `at.risk` for the at-risk functions of `fit.types`.

- verbose.exponential:

  Merged into `verbose`.

- fit.treatment:

  See `treatment_model` in
  [`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md).

- cut.time, cut.one.way, cut.Tk, cut.Tk.values, two.way, reduce.NK,
  hal.sl, screen.two.way, lambda.cvs, event.dependent.cv,
  npenalize.vars, V, seed.hal, reduce.seed.dependence, penalize.time,
  verbose.hal:

  HAL settings, see `hal_control` in
  [`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md).

- use.cores:

  See `cores` in
  [`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md).

- use.cores.prediction, max.Tk, cv.glmnet, browse, browse.hal:

  Ignored.

- return.parameters.for.simulation:

  No longer supported, since `sim.from.data()` has been removed. Use
  [`simevent::sim_model_from_fits()`](https://rdrr.io/pkg/simevent/man/sim_model_from_fits.html)
  instead.

## Value

See
[`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md).
