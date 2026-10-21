# Calibrate a shape-parameter intervention (deprecated)

`calibration.fun()` is deprecated as of MultiStateTMLE 0.1.0. Use
[`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md)
instead.

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
  use.cores = 1,
  target = "outcome",
  z.name = "z",
  ...
)
```

## Arguments

- initial.fit:

  See `initial_fit` in
  [`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md).

- output.eic:

  Ignored;
  [`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md)
  always returns the influence curves.

- tau, tau.z, browse:

  Ignored; the horizon is that of `initial.fit`.

- use.cores:

  See `cores` in
  [`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md).

- target, a, theta, rho, delta, verbose:

  See
  [`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md).
  `theta` is ignored if `rho`, `delta` or `rho.1a` is given.

- z.name, rho.1a:

  See the snake_case equivalents in
  [`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md).

- ...:

  Passed on to
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md),
  with dotted argument names translated to their snake_case equivalents.

## Value

See
[`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md).
