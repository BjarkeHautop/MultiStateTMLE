# TMLE under a shape-parameter intervention (deprecated)

`tmle.alpha.fun()` is deprecated as of MultiStateTMLE 0.1.0. Use
[`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md)
instead.

## Usage

``` r
tmle.alpha.fun(
  target = "z",
  tau = 1.2,
  alpha = 1,
  z.name = "z",
  alpha.list = NULL,
  a = NULL,
  initial.fit = NULL,
  dt = NULL,
  years.lost = NULL,
  only.first = NULL,
  target.only.in.state = NULL,
  target.by.state = FALSE,
  conv.const = 1,
  one.step = FALSE,
  verbose = FALSE,
  max.iter = 10,
  min.iter = 0,
  use.cores = 1,
  truncate.weights = 0,
  output.convergence = FALSE,
  output.eic = FALSE,
  output.weights = NULL,
  output.a.weights = NULL,
  browse = FALSE,
  verbose.exponential = FALSE,
  ...
)
```

## Arguments

- target, alpha, a, verbose:

  See
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md).

- tau:

  Only used when `initial.fit` is `NULL`; the horizon is that of
  `initial.fit`.

- z.name:

  See `z_name` in
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md).
  If `initial.fit` has no such process, `z.name` is ignored, as is
  `alpha` with a warning.

- alpha.list, target.only.in.state, target.by.state, conv.const,
  one.step, max.iter, min.iter, truncate.weights:

  See the snake_case equivalents in
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md)
  and its `control`; `truncate.weights = 0` is no truncation.

- initial.fit:

  See `initial_fit` in
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md).
  If `NULL`, it is computed by
  `prepare.initial(tau = tau, a = a, verbose = verbose, ...)`.

- dt, browse, verbose.exponential:

  Ignored.

- years.lost:

  If given, `estimand = "years_lost"` with
  `control$block_size = years.lost` in
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md).

- only.first:

  See `only_first` in
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md);
  `NULL` means `FALSE`.

- use.cores:

  See `cores` in
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md).

- output.convergence, output.eic, output.weights, output.a.weights:

  Ignored;
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md)
  always returns the convergence diagnostics, influence curve and weight
  summaries.

- ...:

  Passed on to
  [`prepare.initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/prepare.initial.md)
  when `initial.fit` is `NULL`.

## Value

See
[`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md).
