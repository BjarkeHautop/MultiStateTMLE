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
    hazard_model(z ~ L0, limit = 1),
    death ~ L0 + z,
    censoring ~ L0
  ),
  tau = 1
)
fit
#> Initial fit of 150 individuals up to tau = 1
#> 
#>         name     type          model at_risk
#>       <char>   <char>         <char>  <lgcl>
#> 1:         z one_jump         z ~ L0    TRUE
#> 2:     death terminal death ~ L0 + z   FALSE
#> 3: censoring terminal censoring ~ L0   FALSE

est <- tmle_estimate(fit, risk(death, z), scale_intensity(z = 0.5))
est
#> TMLE estimate of risk(death, z)
#> under scale_intensity(z = 0.5)
#> 
#>    target   estimate         se      lower     upper
#>    <char>      <num>      <num>      <num>     <num>
#> 1:  death 0.30405495 0.03853505 0.22852763 0.3795823
#> 2:      z 0.08486244 0.01787691 0.04982433 0.1199005
#> 
#> death: Initial estimate: 0.3124; converged after 2 iterations
#> z: Initial estimate: 0.08486; converged after 1 iterations

# The estimates, named by target process:
coef(est)
#>      death          z 
#> 0.30405495 0.08486244 

# Confidence intervals at another level, for all or some targets:
confint(est, level = 0.9)
#>              5 %      95 %
#> death 0.24067043 0.3674395
#> z     0.05545753 0.1142673
confint(est, "death")
#>           2.5 %    97.5 %
#> death 0.2285276 0.3795823
```
