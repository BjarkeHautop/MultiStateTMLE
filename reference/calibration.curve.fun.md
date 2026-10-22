# Calibration curve over a grid of interventions (deprecated)

`calibration.curve.fun()` is deprecated as of MultiStateTMLE 0.1.0. Use
[`calibration_curve()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibration_curve.md)
instead.

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

  See `initial_fit` in
  [`calibration_curve()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibration_curve.md).

- output.eic:

  Ignored;
  [`calibration_curve()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibration_curve.md)
  always returns the influence curves.

- tau, tau.z, browse:

  Ignored; the horizon is that of `initial.fit`.

- use.cores:

  See `cores` in
  [`calibration_curve()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibration_curve.md).

- target, a, verbose:

  See
  [`calibration_curve()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibration_curve.md).

- z.name, alpha.grid:

  See the snake_case equivalents in
  [`calibration_curve()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibration_curve.md).

- min.iter, target.by.state:

  See the snake_case equivalents in the `control` of
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md).

- block.size.z, block.size.target:

  If given, `exposure_estimand = "years_lost"` and
  `estimand = "years_lost"`, respectively, in
  [`calibration_curve()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibration_curve.md),
  with this `control$block_size`. If both are given, they must be equal.

- ...:

  Passed on to
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md),
  with dotted argument names translated to their snake_case equivalents.

## Value

See
[`calibration_curve()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibration_curve.md).
