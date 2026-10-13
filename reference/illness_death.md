# Simulated illness-death data

Simulated event history data for 150 individuals at risk of an illness
`z`, which can occur at most once, and of death `outcome1`. Both have
exponential hazards that depend on a standard normal baseline covariate
`L0`, the illness raises the hazard of death, and follow-up is subject
to independent censoring. Simulated with
[`simevent::sim_events()`](https://rdrr.io/pkg/simevent/man/sim_events.html),
see `data-raw/illness_death.R` in the package source.

## Usage

``` r
illness_death
```

## Format

A data frame with one row per event and columns:

- id:

  Individual identifier.

- time:

  Event time.

- delta:

  Event type: `0` censoring, `1` death (`outcome1`), `2` illness (`z`).

- L0:

  Baseline covariate.
