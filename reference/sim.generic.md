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
  fallback when those arguments are not supplied.

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
