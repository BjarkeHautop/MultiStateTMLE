# Simulate event history data from parameters fitted to observed data

Simulates new multistate event history data using process/baseline
parameters previously estimated from observed data (e.g. via
[`simevent::simEventCox()`](https://rdrr.io/pkg/simevent/man/simEventCox.html)),
via
[`simevent::simEventData()`](https://rdrr.io/pkg/simevent/man/simEventData.html).

## Usage

``` r
sim.from.data(
  n = 500,
  sim.parameters,
  cens = 1,
  alpha.intervention = list(),
  baseline.intervention = list(),
  browse = FALSE,
  verbose = FALSE
)
```

## Arguments

- n:

  number of individuals to simulate.

- sim.parameters:

  list of fitted simulation parameters, with one entry per process
  (Weibull `weibull.parameters` and `cox.parameters`), a
  `baseline.summary` entry describing baseline covariate distributions,
  and a `model.structure` entry describing process names/types/order.

- cens:

  numeric at-risk indicator scaling for the censoring process.

- alpha.intervention:

  named list of multiplicative interventions on process intensities
  (`eta`), keyed by process name.

- baseline.intervention:

  named list of interventions that fix a baseline covariate to a
  constant value.

- browse:

  logical; if `TRUE`, drop into
  [`browser()`](https://rdrr.io/r/base/browser.html) before simulating.

- verbose:

  logical; if `TRUE`, print the resolved `eta`/`nu`/`beta` parameters.

## Value

A `data.table` of simulated event history data with columns `id`,
`time`, `delta`, baseline covariates, and one column per process.
