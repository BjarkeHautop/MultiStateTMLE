# Compute a calibration curve over a grid of shape-parameter interventions

Evaluates
[`tmle.alpha.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle.alpha.fun.md)
on the `z` process across a grid of shape interventions `alpha.grid`,
and the resulting effect on `target` for each, tracing out a calibration
curve relating `alpha` to the target estimand.

## Usage

``` r
calibration.curve.fun(
  initial.fit = NULL,
  a = NULL,
  alpha.grid = seq(0, 5, length = 10),
  verbose = TRUE,
  output.eic = FALSE,
  tau = 1.2,
  tau.z = tau,
  use.cores = 1,
  target = "outcome",
  z.name = "z",
  min.iter = 1,
  target.by.state = FALSE,
  block.size.z = NULL,
  block.size.target = NULL,
  ...
)
```

## Arguments

- initial.fit:

  initial fit object, as returned by
  [`prepare.initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/prepare.initial.md).

- a:

  optional intervention on the baseline treatment `A0`, passed on to
  [`tmle.alpha.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle.alpha.fun.md).

- alpha.grid:

  numeric vector of shape-parameter interventions to evaluate.

- verbose:

  logical; if `TRUE`, print progress.

- output.eic:

  logical; if `TRUE`, include influence-function values in the output.

- tau:

  follow-up horizon for the `target` process.

- tau.z:

  follow-up horizon for the `z` process (defaults to `tau`).

- use.cores:

  number of cores to use.

- target:

  character; name of the outcome process the curve is computed for.

- z.name:

  character; name of the process on which the shape intervention `alpha`
  acts.

- min.iter:

  minimum number of TMLE update iterations.

- target.by.state:

  logical; whether to compute clever covariates by state, passed on to
  [`tmle.alpha.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle.alpha.fun.md).

- block.size.z:

  optional block size used when evaluating the `z` process.

- block.size.target:

  optional block size used when evaluating `target`.

- ...:

  additional arguments passed on to
  [`tmle.alpha.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle.alpha.fun.md).

## Value

A `data.table` with one row per value in `alpha.grid`, giving the
corresponding TMLE estimate of the effect on `target`.

## Examples

``` r
set.seed(1405)
baseline <- list(L0 = function(N) rnorm(N))
processes <- list(
  z = list(type = "one.jump", eta = 0.2, nu = 1),
  outcome1 = list(type = "terminal", eta = 0.3, nu = 1),
  censoring = list(type = "censoring", eta = 0.1, nu = 1)
)
effects <- list(
  c("L0", "z", 0.5),
  c("z", "outcome1", 0.7)
)
dt <- sim.generic(baseline = baseline, processes = processes, effects = effects, n = 30)

initial.fit <- prepare.initial(
  dt,
  tau = 1,
  fit.types = list(
    z = list(
      model = "Surv(tstart, tstop, delta == 2)~L0",
      fit = "cox",
      at.risk = function(dt) (dt[["z"]] == 0)
    ),
    outcome1 = list(model = "Surv(tstart, tstop, delta == 1)~L0+z", fit = "cox"),
    censoring = list(model = "Surv(tstart, tstop, delta == 0)~L0", fit = "cox")
  ),
  verbose = FALSE
)

# Trace the effect on outcome1 across a small grid of shape interventions.
curve <- calibration.curve.fun(
  initial.fit = initial.fit,
  alpha.grid = c(0, 1),
  tau = 1,
  target = "outcome1",
  z.name = "z",
  use.cores = 1,
  verbose = FALSE
)
curve$estimate
#>    alpha exposure.est exposure.se target.est  target.se
#>    <num>        <num>       <num>      <num>      <num>
#> 1:     0    0.0000000  0.00000000  0.3083278 0.09156753
#> 2:     1    0.2016211  0.09379827  0.3308559 0.09389534
```
