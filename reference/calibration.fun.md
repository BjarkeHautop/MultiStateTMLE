# Calibrate a shape-parameter intervention to a target effect size

Finds the intensity-shape intervention `alpha` (via
[`estimate.alpha.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/estimate.alpha.fun.md)
and
[`tmle.alpha.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle.alpha.fun.md))
for the `z` process such that the resulting effect on `target` matches a
prespecified fraction `rho` of the effect under `alpha = 1`, or a fixed
`theta`/`delta`.

## Usage

``` r
calibration.fun(
  initial.fit = NULL,
  a = NULL,
  theta = 0.5,
  rho = NULL,
  delta = NULL,
  rho.1a = NULL,
  verbose = TRUE,
  output.eic = FALSE,
  tau = 1.2,
  tau.z = tau,
  use.cores = 1,
  target = "outcome",
  z.name = "z",
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

- theta:

  optional numeric; target value for the estimand to calibrate to
  directly.

- rho:

  optional numeric; target fraction of the `alpha = 1` effect on
  `target` to calibrate to.

- delta:

  optional numeric; target absolute difference from the `alpha = 1`
  effect on `target` to calibrate to.

- rho.1a:

  optional numeric; as `rho`, but relative to the effect under `a` fixed
  to 1 rather than `alpha = 1`.

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

  character; name of the outcome process being calibrated to.

- z.name:

  character; name of the process on which the shape intervention `alpha`
  acts.

- ...:

  additional arguments passed on to
  [`tmle.alpha.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle.alpha.fun.md).

## Value

A list with the calibrated `alpha` (as found by
[`estimate.alpha.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/estimate.alpha.fun.md))
and the corresponding TMLE fit(s).

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

# Calibrate alpha so the risk of outcome1 by tau = 1 equals 0.2.
cal <- calibration.fun(
  initial.fit = initial.fit,
  theta = 0.2,
  tau = 1,
  target = "outcome1",
  z.name = "z",
  use.cores = 1,
  verbose = FALSE
)
cal$estimate
#>       alpha.est        alpha.se      target.est       target.se target.se.crude 
#>      1.00000000      0.54173579      0.34707334      0.09416282      0.09024551 
#>           theta        theta.se 
#>      0.20000000      0.00000000 
```
