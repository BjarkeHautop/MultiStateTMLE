# Specify the hazard model of an event process

Describes how
[`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md)
models the hazard of one event process. A plain formula such as
`death ~ L0 + z` in `models` is short for
`hazard_model(death ~ L0 + z)`; use `hazard_model()` to also restrict
when an individual is at risk of the process, or to change how it is
fitted.

## Usage

``` r
hazard_model(formula, limit = NULL, at_risk = NULL, fit = c("cox", "hal"))
```

## Arguments

- formula:

  A two-sided formula. The left-hand side names the process, matching
  its name in the `event` column of the data, and the right-hand side
  gives the covariates of its hazard, as in a
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)
  formula. Other processes can be used as covariates by name, giving
  their number of events so far.

- limit:

  Optional positive integer. Maximum number of events of the process per
  individual: individuals are only at risk while they have had fewer
  than `limit` events. `limit = 1` is an event that can happen at most
  once, as `sim_process(limit = 1)` in simevent.

- at_risk:

  Optional one-sided formula, such as `~ z == 0`, of a condition on the
  processes (their number of events so far). Individuals are only at
  risk of the process while it holds. Combined with `limit` if both are
  given.

- fit:

  How the hazard is fitted: `"cox"` (default), a Cox model, or `"hal"`,
  the highly adaptive lasso (experimental, see
  [`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md)).

## Value

An object of class `hazard_model`.

## Examples

``` r
# A plain formula in fit_initial()'s `models` is short for
# hazard_model(formula):
hazard_model(death ~ L0 + z)
#> Hazard model: death ~ L0 + z

# An illness that can happen at most once:
hazard_model(z ~ L0, limit = 1)
#> Hazard model: z ~ L0
#>   limit: 1

# A relapse that can happen at most twice, and only before death of
# another cause, w:
hazard_model(relapse ~ L0, limit = 2, at_risk = ~ w == 0)
#> Hazard model: relapse ~ L0
#>   limit: 2
#>   at risk while: w == 0

# Death depending on having had two or more relapses. As in any R model
# formula, wrap arithmetic and comparisons in I():
hazard_model(death ~ L0 + I(relapse >= 2))
#> Hazard model: death ~ L0 + I(relapse >= 2)

# Use them in fit_initial(), mixed with plain formulas:
fit <- fit_initial(
  illness_death,
  models = list(
    hazard_model(z ~ L0, limit = 1),
    death ~ L0 + z,
    censoring ~ L0
  ),
  tau = 1
)
fit$processes
#>         name     type          model at_risk
#>       <char>   <char>         <char>  <lgcl>
#> 1:         z one_jump         z ~ L0    TRUE
#> 2:     death terminal death ~ L0 + z   FALSE
#> 3: censoring terminal censoring ~ L0   FALSE
```
