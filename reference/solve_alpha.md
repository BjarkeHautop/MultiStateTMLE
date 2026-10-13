# Solve for the intervention achieving a target value

Finds the intensity intervention `alpha` at which the estimate returned
by `fun(alpha)` equals `theta`, assuming the estimate is increasing in
`alpha`. The search first expands a bracket around the solution, then
narrows it by secant steps (with bisection as a fallback), caching
evaluations of `fun` along the way.

## Usage

``` r
solve_alpha(
  fun,
  theta,
  c_n,
  alpha_init = 1,
  expand_up = 2,
  expand_down = 0.25,
  max_iter = 100,
  alpha_min = 0.001,
  alpha_max = 100,
  verbose = FALSE,
  trace_every = 1L
)
```

## Arguments

- fun:

  A function of `alpha` returning a list with elements `estimate` (a
  named vector containing `tmle.est` or `one.step.est`, and `se`) and
  `eic`, such as a call to
  [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md).

- theta:

  Number. The target value of the estimate.

- c_n:

  Positive number. Tolerance; the search stops once the estimate is
  within `c_n` of `theta`.

- alpha_init:

  Number. Starting value of `alpha`.

- expand_up:

  Number greater than 1. Factor by which `alpha` is multiplied while
  searching upwards for a bracket.

- expand_down:

  Number between 0 and 1. Factor by which `alpha` is multiplied while
  searching downwards for a bracket.

- max_iter:

  Positive integer. Maximum number of steps in each of the bracketing
  and narrowing phases.

- alpha_min, alpha_max:

  Numbers with `0 < alpha_min < alpha_max`. Bounds on `alpha`.

- verbose:

  Logical. If `TRUE`, print progress every `trace_every` evaluations.

- trace_every:

  Positive integer. How often to print progress when `verbose = TRUE`.

## Value

A list with elements `alpha.hat` (the best `alpha` found), `converged`
(whether it is within `c_n` of `theta`), `dist` (its distance to
`theta`), `c.n`, `grid` (a data frame of all evaluated `alpha`, `psi`
and `se`) and `eic` (the influence curve at `alpha.hat`).

## Details

This solves the calibration equation of Rytgaard and van der Laan, and
is used by
[`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md).

## Examples

``` r
# psi(alpha) = 0.6 * (1 - exp(-alpha)) has a known inverse, so
# solve_alpha() should recover alpha = -log(1 - 0.3 / 0.6).
fun <- function(alpha) {
  psi <- 0.6 * (1 - exp(-alpha))
  list(estimate = c(tmle.est = psi, se = 0.01), eic = rep(0, 10))
}
res <- solve_alpha(fun, theta = 0.3, c_n = 1e-6)
res$alpha.hat
#> [1] 0.6931482
-log(1 - 0.3 / 0.6)
#> [1] 0.6931472
```
