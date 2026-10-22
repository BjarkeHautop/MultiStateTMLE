# Methods for the result objects

Print, extract the estimates of, and compute confidence intervals for
the objects returned by
[`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md),
[`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md),
[`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md),
[`calibration_curve()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibration_curve.md)
and
[`calibrated_contrasts()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrated_contrasts.md).

## Usage

``` r
# S3 method for class 'initial_fit'
print(x, ...)

# S3 method for class 'tmle_alpha'
print(x, ...)

# S3 method for class 'calibrate_alpha'
print(x, ...)

# S3 method for class 'calibration_curve'
print(x, ...)

# S3 method for class 'calibrated_contrasts'
print(x, ...)

# S3 method for class 'tmle_alpha'
coef(object, ...)

# S3 method for class 'calibrate_alpha'
coef(object, ...)

# S3 method for class 'tmle_alpha'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'calibrate_alpha'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- x, object:

  A result object.

- ...:

  Ignored.

- parm:

  Optional character vector of targets to compute confidence intervals
  for. Defaults to all.

- level:

  Number between 0 and 1. Confidence level.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.
[`coef()`](https://rdrr.io/r/stats/coef.html) returns a numeric vector
of estimates named by target, and
[`confint()`](https://rdrr.io/r/stats/confint.html) a matrix of Wald
confidence intervals with one row per target.

## Examples

``` r
fit <- fit_initial(
  illness_death,
  fit_types = list(
    z = list(
      delta = 2,
      model = ~L0,
      at_risk = function(dt) dt[["z"]] == 0
    ),
    death = list(delta = 1, model = ~ L0 + z),
    censoring = list(delta = 0, model = ~L0)
  ),
  tau = 1
)
fit
#> Initial fit of 150 individuals up to tau = 1
#> 
#>         name delta     type   model at_risk
#>       <char> <num>   <char>  <char>  <lgcl>
#> 1:         z     2 one_jump     ~L0    TRUE
#> 2:     death     1 terminal ~L0 + z   FALSE
#> 3: censoring     0 terminal     ~L0   FALSE

est <- tmle_alpha(fit, target = "death", alpha = 0.5, z_name = "z")
est
#> TMLE estimate under alpha = 0.5 on z
#> 
#>    target  estimate         se     lower     upper
#>    <char>     <num>      <num>     <num>     <num>
#> 1:  death 0.3040549 0.03853505 0.2285276 0.3795823
#> 
#> Initial estimate: 0.3124; converged after 2 iterations
coef(est)
#>     death 
#> 0.3040549 
confint(est, level = 0.9)
#>             5 %      95 %
#> death 0.2406704 0.3674395
```
