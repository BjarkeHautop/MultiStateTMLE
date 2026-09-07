# Estimate the numerical derivative of an alpha-indexed estimand

Approximates the derivative with respect to `alpha` of the TMLE/one-step
estimate returned by `fun`, using a symmetric finite-difference at
`alpha_hat`.

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

  numeric; the point at which to estimate the derivative.

- fun:

  function taking `alpha`, `parameter`, and `...`, returning a list with
  an `estimate` element containing `tmle.est`/`one.step.est`.

- parameter:

  character; the target process/outcome passed to `fun`.

- h:

  numeric step size for the finite difference; shrunk if it exceeds
  `alpha_hat`.

- ...:

  additional arguments passed on to `fun`.

## Value

numeric; the estimated derivative.
