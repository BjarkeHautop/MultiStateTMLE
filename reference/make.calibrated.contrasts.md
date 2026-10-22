# Total, mediated and calibrated contrasts (deprecated)

`make.calibrated.contrasts()` is deprecated as of MultiStateTMLE 0.1.0.
Use
[`calibrated_contrasts()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrated_contrasts.md)
instead.

## Usage

``` r
make.calibrated.contrasts(
  target = NULL,
  treatment.fit,
  placebo.fit,
  calibrated.fit,
  conf.level = 0.95
)
```

## Arguments

- target:

  If given, must be the target of `treatment.fit`, which
  [`calibrated_contrasts()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrated_contrasts.md)
  always uses.

- treatment.fit, placebo.fit, calibrated.fit:

  See the snake_case equivalents in
  [`calibrated_contrasts()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrated_contrasts.md).

- conf.level:

  See `level` in
  [`calibrated_contrasts()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrated_contrasts.md).

## Value

See
[`calibrated_contrasts()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrated_contrasts.md).
