# Methods for the result objects

Print, extract the estimates of, and compute confidence intervals for
the objects returned by
[`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md),
[`tmle_estimate()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_estimate.md),
[`intervention_curve()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/intervention_curve.md)
and
[`calibrated_contrasts()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrated_contrasts.md).

## Usage

``` r
# S3 method for class 'initial_fit'
print(x, ...)

# S3 method for class 'tmle_estimate'
print(x, ...)

# S3 method for class 'intervention_curve'
print(x, ...)

# S3 method for class 'calibrated_contrasts'
print(x, ...)

# S3 method for class 'tmle_estimate'
coef(object, ...)

# S3 method for class 'tmle_estimate'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- x, object:

  A result object.

- ...:

  Ignored.

- parm:

  Optional character vector of target processes to compute confidence
  intervals for. Defaults to all.

- level:

  Number between 0 and 1. Confidence level.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.
[`coef()`](https://rdrr.io/r/stats/coef.html) returns a numeric vector
of estimates named by target process, and
[`confint()`](https://rdrr.io/r/stats/confint.html) a matrix of Wald
confidence intervals with one row per target process.

## Examples

``` r
fit <- fit_initial(
  illness_death,
  models = list(
    hazard_model(illness ~ biomarker + treatment, limit = 1),
    death ~ biomarker + treatment + illness,
    censoring ~ biomarker
  ),
  tau = 1
)
fit
#> Initial fit of 200 individuals up to tau = 1
#> 
#>         name     type                                   model at_risk tracked
#>       <char>   <char>                                  <char>  <lgcl>   <int>
#> 1:   illness one_jump         illness ~ biomarker + treatment    TRUE       1
#> 2:     death terminal death ~ biomarker + treatment + illness   FALSE      NA
#> 3: censoring terminal                   censoring ~ biomarker   FALSE      NA

est <- tmle_estimate(
  fit,
  risk(death, illness),
  scale_intensity(illness = 0.5)
)
est
#> TMLE estimate of risk(death, illness)
#> under scale_intensity(illness = 0.5)
#> 
#>     target  estimate         se     lower     upper
#>     <char>     <num>      <num>     <num>     <num>
#> 1:   death 0.2906665 0.03341164 0.2251809 0.3561521
#> 2: illness 0.2169103 0.02428345 0.1693156 0.2645050
#> 
#> death: Initial estimate: 0.2907; no targeting needed, as it already solves the EIC equation
#> illness: Initial estimate: 0.2169; no targeting needed, as it already solves the EIC equation

# The estimates, named by target process:
coef(est)
#>     death   illness 
#> 0.2906665 0.2169103 

# Confidence intervals at another level, for all or some targets:
confint(est, level = 0.9)
#>               5 %      95 %
#> death   0.2357092 0.3456237
#> illness 0.1769676 0.2568531
confint(est, "death")
#>           2.5 %    97.5 %
#> death 0.2251809 0.3561521
```
