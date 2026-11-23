# Build a simulation model from an initial fit

Turns the Cox models of an
[`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md)
fit into a
[`simevent::sim_model()`](https://rdrr.io/pkg/simevent/man/sim_model.html),
with
[`simevent::sim_model_from_fits()`](https://rdrr.io/pkg/simevent/man/sim_model_from_fits.html),
so data resembling the observed data can be simulated with
[`simevent::simulate_events()`](https://rdrr.io/pkg/simevent/man/simulate_events.html),
e.g. to check an estimator where the truth is known. The process types
and limits are taken from the fit.

## Usage

``` r
as_sim_model(fit, data)
```

## Arguments

- fit:

  An initial fit, as returned by
  [`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md).

- data:

  The `data` that `fit` was fitted to, to resample covariates from.

## Value

A
[`simevent::sim_model()`](https://rdrr.io/pkg/simevent/man/sim_model.html).

## Details

Covariates are resampled from `data`, one row per individual. A process
that a hazard depends on enters it as in the fit: a recurrent process as
the indicators of at least 1, ..., `max_count` jumps, and any other as
its number of jumps so far. A process that jumps at most once in the
data (type `"one_jump"`) can jump at most once in the simulation too.

Fits with an `at_risk` condition in a
[`hazard_model()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/hazard_model.md),
or with `event_time()` covariates, are not supported.

## Examples

``` r
# \donttest{
if (requireNamespace("simevent", quietly = TRUE)) {
  fit <- fit_initial(
    illness_death,
    models = list(
      hazard_model(illness ~ biomarker + treatment, limit = 1),
      death ~ biomarker + treatment + illness,
      censoring ~ biomarker
    ),
    tau = 1
  )
  model <- as_sim_model(fit, illness_death)
  simulated <- simevent::simulate_events(model, n = 500, seed = 1)
  head(simulated)
}
#> Key: <id>
#>       id      time   event   biomarker treatment illness
#>    <int>     <num>  <fctr>       <num>     <num>   <num>
#> 1:     1 0.1730469 illness  0.03797316         0       1
#> 2:     1 0.7092549   death  0.03797316         0       1
#> 3:     2 2.5093113   death  2.94325166         1       0
#> 4:     3 0.8452039 illness -0.17643927         1       1
#> 5:     3 0.9287521   death -0.17643927         1       1
#> 6:     4 1.2212702   death -2.16792215         0       0
# }
```
