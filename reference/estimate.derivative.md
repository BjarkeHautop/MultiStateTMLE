# Estimate the derivative with respect to alpha (deprecated)

`estimate.derivative()` is deprecated as of MultiStateTMLE 0.1.0. Use
[`estimate_derivative()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/estimate_derivative.md)
instead.

## Usage

``` r
estimate.derivative(
  alpha_hat,
  fun = alpha.est.fun,
  parameter = "z",
  h = length(unique(dt[, id]))^{
-1/6
 },
  ...
)
```

## Arguments

- alpha_hat:

  See `alpha` in
  [`estimate_derivative()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/estimate_derivative.md).

- fun, h:

  See
  [`estimate_derivative()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/estimate_derivative.md).

- parameter:

  Passed on to `fun` as `parameter`.

- ...:

  Passed on to `fun`.

## Value

See
[`estimate_derivative()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/estimate_derivative.md).
