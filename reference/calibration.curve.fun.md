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
  browse = FALSE,
  verbose = TRUE,
  output.eic = FALSE,
  tau = 1.2,
  tau.z = tau,
  use.cores = 50,
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

- browse:

  logical; if `TRUE`, drop into
  [`browser()`](https://rdrr.io/r/base/browser.html).

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
