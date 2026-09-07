# Targeted maximum likelihood estimation under a shape-parameter intervention

Estimates the effect of a stochastic intensity intervention (scaling the
shape of the `z` process by `alpha`, optionally combined with an
intervention `a` on baseline treatment) on `target`, via an iterative
TMLE update of the initial fit produced by
[`prepare.initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/prepare.initial.md).

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
  use.cores = 50,
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

- target:

  character; name of the target process/outcome.

- tau:

  follow-up horizon (in time units) for `target`.

- alpha:

  numeric; multiplicative shape-parameter intervention applied to the
  `z` process (ignored if `alpha.list` is supplied).

- z.name:

  character; name of the process on which `alpha` acts.

- alpha.list:

  optional named list of intervention functions, one per process, of the
  form `function(time, covariates)`; supersedes `alpha`.

- a:

  optional intervention on the baseline treatment `A0` (must match the
  `a` used to build `initial.fit`, if `initial.fit` is supplied).

- initial.fit:

  initial fit object, as returned by
  [`prepare.initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/prepare.initial.md);
  computed from `dt`/`tau`/`a`/`...` if not supplied.

- dt:

  data.table of observed event history data, passed to
  [`prepare.initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/prepare.initial.md)
  when `initial.fit` is not supplied.

- years.lost:

  optional; if supplied, computes years-lost-type estimands with this
  block size.

- only.first:

  logical; if `TRUE`, restrict `target` to its first occurrence (treat
  it as a one-jump process).

- target.only.in.state:

  optional function of `states` identifying the subset of states in
  which `target` events are counted.

- target.by.state:

  logical; whether to compute clever covariates by state rather than
  pooled.

- conv.const:

  constant scaling the convergence criterion on the mean efficient
  influence curve.

- one.step:

  logical; if `TRUE`, return a one-step estimator instead of iterating
  to TMLE convergence.

- verbose:

  logical; if `TRUE`, print iteration progress.

- max.iter:

  maximum number of TMLE update iterations.

- min.iter:

  minimum number of TMLE update iterations before checking convergence.

- use.cores:

  number of cores to use for the per-id computations.

- truncate.weights:

  numeric; if `> 0`, truncate clever weights above this value.

- output.convergence:

  logical; if `TRUE`, include convergence diagnostics in the output.

- output.eic:

  logical; if `TRUE`, include the efficient influence curve in the
  output.

- output.weights:

  optional numeric vector of quantiles at which to summarize the
  clever/censoring weights.

- output.a.weights:

  optional numeric vector of quantiles at which to summarize the
  `a`-specific clever weights.

- browse:

  logical; if `TRUE`, drop into
  [`browser()`](https://rdrr.io/r/base/browser.html).

- verbose.exponential:

  logical; if `TRUE`, print diagnostics from the exponential TMLE update
  step.

- ...:

  additional arguments passed on to
  [`prepare.initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/prepare.initial.md)
  when `initial.fit` is not supplied.

## Value

A list with the TMLE `estimate`, standard error, and (depending on the
arguments above) the efficient influence curve, convergence diagnostics,
and weight summaries.
