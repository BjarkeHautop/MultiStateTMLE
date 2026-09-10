# Search for the alpha achieving a target estimand value

Searches over the shape-parameter intervention `alpha` for the value at
which `fun(alpha, ...)`'s estimate equals `theta`, via bracket expansion
followed by root finding, caching evaluations of `fun` along the way.

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
  use.cores = 50,
  verbose = FALSE,
  trace_every = 1L
)
```

## Arguments

- theta:

  numeric; the target value of the estimand to solve for.

- fun:

  function of `alpha` (plus `...`) returning a list with an `estimate`
  element and an `eic` element, e.g.
  [`tmle.alpha.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle.alpha.fun.md).

- c_n:

  numeric; a finite tolerance/scale constant used to judge convergence
  (required, no default).

- alpha_init:

  numeric; initial value of `alpha` to evaluate.

- expand_up:

  multiplicative factor used to expand the bracket upward.

- expand_down:

  multiplicative factor used to expand the bracket downward.

- max_iter:

  maximum number of root-finding iterations.

- alpha_min:

  lower bound for `alpha` (must be `> 0`).

- alpha_max:

  upper bound for `alpha` (must be `> alpha_min`).

- use.cores:

  number of cores to pass on to `fun`.

- verbose:

  logical; if `TRUE`, print progress at every `trace_every` evaluations.

- trace_every:

  integer; frequency (in evaluations) of progress output when
  `verbose = TRUE`.

## Value

A list describing the solution, including the selected `alpha` and the
corresponding evaluation of `fun`.

## Examples

``` r
# Psi_z^alpha(P) = 0.6 * (1 - exp(-alpha)), with a known inverse,
# used here to check that estimate.alpha.fun recovers it numerically.
# eic mimics a real per-subject efficient influence curve, as returned
# by tmle.alpha.fun().
set.seed(1405)
fun <- function(alpha) {
  psi <- 0.6 * (1 - exp(-alpha))
  n <- 200
  ic <- rnorm(n, sd = 0.05)
  ic <- ic - mean(ic)
  list(estimate = c(tmle.est = psi, se = sd(ic) / sqrt(n)), eic = ic)
}

res <- estimate.alpha.fun(theta = 0.3, fun = fun, c_n = 1e-4, alpha_init = 1)
res$alpha.hat
#> [1] 0.6932652
res$converged
#> [1] TRUE
```
