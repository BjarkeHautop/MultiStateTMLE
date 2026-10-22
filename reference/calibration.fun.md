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

- theta, rho, delta, rho.1a:

  The calibration target: `theta` with `scale = "absolute"`, `rho` with
  `scale = "ratio"`, `delta` with `scale = "difference"`, or `rho.1a`
  with `scale = "ratio"` and `reference_a = 1 - a` in
  [`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md).
  `theta` is ignored if one of the others is given.

- output.eic:

  Ignored;
  [`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md)
  always returns the influence curves.

- tau, tau.z, browse:

  Ignored; the horizon is that of `initial.fit`.

- use.cores:

  See `cores` in
  [`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md).

- target, a, verbose:

  See
  [`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md).

- z.name:

  See `z_name` in
  [`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md).

- ...:

  Passed on to
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md),
  with dotted argument names translated to their snake_case equivalents.

## Value

See
[`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md).
