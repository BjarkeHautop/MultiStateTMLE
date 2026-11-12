# Target parameters

Describe what
[`tmle_estimate()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_estimate.md)
and
[`intervention_curve()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/intervention_curve.md)
estimate for one or more event processes, by the horizon `tau` of the
initial fit:

- `risk()`: the risk of an event by `tau`, \\E\[N(\tau)\]\\, which is
  the expected number of events for a recurrent process.

- `years_lost()`: the expected time lost before `tau`, the integral of
  the risk up to `tau`.

## Usage

``` r
risk(..., when = NULL, only_first = FALSE)

years_lost(..., when = NULL, only_first = FALSE)
```

## Arguments

- ...:

  The event processes, as bare names (`death`) or strings (`"death"`).
  Each is estimated separately, under the same intervention. To give the
  names in a character vector `x`, use `do.call(risk, as.list(x))`.

- when:

  Optional one-sided formula of a condition on the state, such as
  `~ z == 0`: only events while it holds count. It may use the processes
  that hazards depend on (their number of events so far).

- only_first:

  Logical. If `TRUE`, only the first event of each process counts.

## Value

An object of class `tmle_target`.

## Examples

``` r
# The risk of death by tau, and the expected years lost before it:
risk(death)
#> Target: risk(death)
years_lost(death)
#> Target: years_lost(death)

# Several processes, each estimated separately under the same intervention:
risk(death, z)
#> Target: risk(death, z)

# Only deaths before illness count:
risk(death, when = ~ z == 0)
#> Target: risk(death, when = ~z == 0)

# For a recurrent process, the risk of a first event rather than the
# expected number of events:
risk(relapse, only_first = TRUE)
#> Target: risk(relapse, only_first = TRUE)

# Processes held in a character vector:
targets <- c("death", "z")
do.call(risk, as.list(targets))
#> Target: risk(death, z)

# Use them in tmle_estimate():
fit <- fit_initial(
  illness_death,
  models = list(
    hazard_model(z ~ L0, limit = 1),
    death ~ L0 + z,
    censoring ~ L0
  ),
  tau = 1
)
tmle_estimate(fit, risk(death, z))
#> TMLE estimate of risk(death, z)
#> without intervention
#> 
#>    target  estimate         se      lower     upper
#>    <char>     <num>      <num>      <num>     <num>
#> 1:  death 0.3250469 0.03906316 0.24848453 0.4016093
#> 2:      z 0.1589128 0.03132140 0.09752402 0.2203016
#> 
#> death: Initial estimate: 0.325; converged after 1 iterations
#> z: Initial estimate: 0.1589; converged after 1 iterations
tmle_estimate(fit, risk(death, when = ~ z == 0))$estimate
#>    target  estimate         se     lower     upper
#>    <char>     <num>      <num>     <num>     <num>
#> 1:  death 0.2608621 0.03686855 0.1886011 0.3331231
```
