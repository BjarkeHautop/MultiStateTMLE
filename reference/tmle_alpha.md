# Targeted estimation under a stochastic intensity intervention

Estimates a target parameter under a stochastic intensity intervention
that scales the intensity of the process `z_name` by `alpha`, optionally
combined with an intervention `a` on the baseline treatment. The initial
fit from
[`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md)
is updated by iterative targeted maximum likelihood estimation (TMLE).

## Usage

``` r
tmle_alpha(
  initial_fit,
  target,
  alpha = 1,
  z_name = "z",
  alpha_list = NULL,
  a = NULL,
  years_lost = NULL,
  only_first = FALSE,
  target_only_in_state = NULL,
  target_by_state = FALSE,
  conv_const = 1,
  one_step = FALSE,
  max_iter = 10,
  min_iter = 0,
  truncate_weights = 0,
  output_convergence = FALSE,
  output_eic = FALSE,
  output_weights = NULL,
  output_a_weights = NULL,
  cores = 1,
  verbose = FALSE
)
```

## Arguments

- initial_fit:

  An initial fit, as returned by
  [`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md).

- target:

  String. Name of the target process.

- alpha:

  Non-negative number. Multiplicative intervention on the intensity of
  `z_name`. `alpha = 1` is no intervention. Ignored if `alpha_list` is
  given.

- z_name:

  String. Name of the process whose intensity `alpha` scales. If no such
  process exists in `initial_fit`, `alpha` is ignored with a warning.

- alpha_list:

  Experimental. Optional named list of functions, one per process,
  giving a time-varying intervention on that process's intensity. Each
  function's arguments name the columns of `initial_fit$tmp.long` it
  depends on. Supersedes `alpha` and `z_name`.

- a:

  Optional baseline treatment value (`0` or `1`) to intervene on. Must
  be one of the values of `a` given to
  [`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md).

- years_lost:

  Optional positive integer. If given, the target parameter is the
  expected time lost (the integral of the risk of `target` up to `tau`)
  instead of the risk at `tau`. The integral is approximated with blocks
  of `years_lost` time points.

- only_first:

  Logical. If `TRUE`, only the first `target` event counts, and only the
  first jump of `z_name` is intervened on.

- target_only_in_state:

  Optional function of the state space (`initial_fit$depend.matrix`)
  returning a logical vector, marking the states in which `target`
  events count.

- target_by_state:

  Logical. If `TRUE`, use state-specific clever covariates in the
  targeting step.

- conv_const:

  Positive number. The TMLE iterations stop when the absolute mean of
  the efficient influence curve falls below `conv_const * se / log(n)`.

- one_step:

  Logical. If `TRUE`, return the one-step estimator instead of iterating
  the TMLE update.

- max_iter:

  Positive integer. Maximum number of TMLE iterations.

- min_iter:

  Non-negative integer. Minimum number of TMLE iterations before
  checking convergence.

- truncate_weights:

  Non-negative number. If positive, the clever weights of the
  intervention are truncated at this value.

- output_convergence:

  Logical. If `TRUE`, include convergence diagnostics in the output.

- output_eic:

  Logical. If `TRUE`, include the efficient influence curve in the
  output.

- output_weights:

  Optional numeric vector of probabilities at which to report quantiles
  of the per-individual maximal intervention and censoring weights.

- output_a_weights:

  Optional numeric vector of probabilities at which to report quantiles
  of the per-individual maximal treatment weights.

- cores:

  Positive integer. Number of cores used by
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html) for
  the per-individual computations.

- verbose:

  Logical. If `TRUE`, print progress.

## Value

A list with element `estimate`, a named vector with the estimate
(`tmle.est`, or `one.step.est` if `one_step = TRUE`), its standard error
`se`, and the initial plug-in estimate `g.est`. Depending on the
arguments, the list also contains `convergence`, `weights`,
`cens.weights`, `a.weights` and `eic`.

## Details

By default the target parameter is the risk of a `target` event by the
horizon `tau` given to
[`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md)
(or the expected number of events for a recurrent `target`).

## Examples

``` r
fit <- fit_initial(
  illness_death,
  fit_types = list(
    z = list(
      model = "Surv(tstart, tstop, delta == 2) ~ L0",
      at_risk = function(dt) dt[["z"]] == 0
    ),
    outcome1 = "Surv(tstart, tstop, delta == 1) ~ L0 + z",
    censoring = "Surv(tstart, tstop, delta == 0) ~ L0"
  ),
  tau = 1
)

# Risk of outcome1 by tau = 1 if the intensity of z were halved.
est <- tmle_alpha(fit, target = "outcome1", alpha = 0.5, z_name = "z")
est$estimate
#>   tmle.est         se      g.est 
#> 0.30405495 0.03853505 0.31239982 
```
