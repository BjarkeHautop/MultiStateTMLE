# Search for the alpha achieving a target value (deprecated)

`estimate.alpha.fun()` is deprecated as of MultiStateTMLE 0.1.0. Use
[`solve_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/solve_alpha.md)
instead.

## Usage

``` r
estimate.alpha.fun(
  theta,
  fun,
  c_n,
  alpha_init = 1,
  expand_up = 2,
  expand_down = 0.25,
  max_iter = 100,
  alpha_min = 0.001,
  alpha_max = 100,
  use.cores = 1,
  verbose = FALSE,
  trace_every = 1L
)
```

## Arguments

- theta, fun, c_n, alpha_init, expand_up, expand_down, max_iter,
  alpha_min, alpha_max, verbose, trace_every:

  See
  [`solve_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/solve_alpha.md).

- use.cores:

  Ignored.

## Value

See
[`solve_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/solve_alpha.md).
