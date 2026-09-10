# Simulate event history data from a generic process specification

Simulates multistate event history data from a set of user-specified
baseline covariates, event processes (with Weibull intensities and
Cox-type effects), and their effects on one another, via
[`simevent::simEventData()`](https://rdrr.io/pkg/simevent/man/simEventData.html).

## Usage

``` r
sim.generic(
  baseline = list(),
  processes = list(),
  effects = list(),
  sim.object = list(),
  cens = 1,
  alpha.intervention = list(),
  baseline.intervention = list(),
  n = 500,
  browse = FALSE
)
```

## Arguments

- baseline:

  named list of baseline covariate generators.

- processes:

  named list of process specifications; each entry has at least a `type`
  (one of `"censoring"`, `"terminal"`, `"one.jump"`, or recurrent), and
  Weibull intensity parameters `eta`/`nu`.

- effects:

  list of `c(from, to, coefficient)` triples specifying Cox effects
  between baseline/process variables and process intensities.

- sim.object:

  optional list with `baseline`/`processes`/`effects` entries, used as a
  fallback when none of `baseline`, `processes`, or `effects` are
  supplied directly.

- cens:

  numeric at-risk indicator scaling for the censoring process.

- alpha.intervention:

  named list of multiplicative interventions on process intensities
  (`eta`), keyed by process name.

- baseline.intervention:

  named list of interventions that fix a baseline covariate to a
  constant value.

- n:

  number of individuals to simulate.

- browse:

  logical; if `TRUE`, drop into
  [`browser()`](https://rdrr.io/r/base/browser.html) before simulating.

## Value

A `data.table` of simulated event history data with columns `id`,
`time`, `delta`, baseline covariates, and one column per process.

## Examples

``` r
# An illness-death process: z is a one-jump "illness" event that raises
# the hazard of the terminal process outcome1, with independent censoring.
set.seed(1405)
baseline <- list(L0 = function(N) rnorm(N))
processes <- list(
  z = list(type = "one.jump", eta = 0.2, nu = 1),
  outcome1 = list(type = "terminal", eta = 0.3, nu = 1),
  censoring = list(type = "censoring", eta = 0.1, nu = 1)
)
effects <- list(
  c("L0", "z", 0.5),
  c("z", "outcome1", 0.7)
)

dt <- sim.generic(baseline = baseline, processes = processes, effects = effects, n = 100)
head(dt)
#> Key: <id>
#>       id      time delta         L0     z
#>    <int>     <num> <int>      <num> <num>
#> 1:     1 1.0836890     2  0.2724785     1
#> 2:     1 2.7886939     1  0.2724785     1
#> 3:     2 0.6539194     1  0.3572619     0
#> 4:     3 0.1212156     1 -0.8616620     0
#> 5:     4 1.1416233     2  0.8083350     1
#> 6:     4 2.6092612     1  0.8083350     1
```
