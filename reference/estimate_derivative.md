# Estimate the derivative of an estimate with respect to alpha

Approximates the derivative with respect to `alpha` of the estimate
returned by `fun`, by a symmetric finite difference at `alpha`.

## Usage

``` r
estimate_derivative(fun, alpha, h, ...)
```

## Arguments

- fun:

  A function of `alpha` (and `...`) returning a list with an element
  `estimate`, a named vector containing `tmle.est` or `one.step.est`,
  such as a call to
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md).

- alpha:

  Positive number. The point at which to estimate the derivative.

- h:

  Positive number. Step size of the finite difference. It is shrunk to
  `3 * alpha / 4` if it exceeds `alpha`, so that `fun` is only evaluated
  at positive values.

- ...:

  Further arguments passed on to `fun`.

## Value

A number, the estimated derivative.

## Examples

``` r
# The derivative of alpha^2 at alpha = 3 is 6.
fun <- function(alpha) list(estimate = c(tmle.est = alpha^2))
estimate_derivative(fun, alpha = 3, h = 0.01)
#> [1] 6
```
