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

- verbose:

  logical; if `TRUE`, print the resolved `eta`/`nu`/`beta` parameters.

## Value

A `data.table` of simulated event history data with columns `id`,
`time`, `delta`, baseline covariates, and one column per process.

## Examples

``` r
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

# Fit Cox/Weibull parameters from observed data, then simulate new data
# from those fitted parameters.
sim.parameters <- prepare.initial(
  dt,
  tau = 1,
  fit.types = list(
    z = list(
      model = "Surv(tstart, tstop, delta == 2)~L0",
      fit = "cox",
      at.risk = function(dt) (dt[["z"]] == 0)
    ),
    outcome1 = list(model = "Surv(tstart, tstop, delta == 1)~L0+z", fit = "cox"),
    censoring = list(model = "Surv(tstart, tstop, delta == 0)~L0", fit = "cox")
  ),
  verbose = FALSE,
  return.parameters.for.simulation = TRUE
)

new.dt <- sim.from.data(n = 100, sim.parameters = sim.parameters)
head(new.dt)
#> Key: <id>
#>       id      time delta        L0     z
#>    <int>     <num> <int>     <num> <num>
#> 1:     1 0.1508955     1 0.2366780     0
#> 2:     2 2.1064530     0 0.8048426     0
#> 3:     3 1.1223962     1 0.9227430     0
#> 4:     4 0.6157989     2 0.1345364     1
#> 5:     4 0.7253412     1 0.1345364     1
#> 6:     5 3.2244288     1 1.6101867     0
```
