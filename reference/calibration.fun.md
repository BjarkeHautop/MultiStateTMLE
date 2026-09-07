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
  browse = FALSE,
  verbose = TRUE,
  output.eic = FALSE,
  tau = 1.2,
  tau.z = tau,
  use.cores = 50,
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
