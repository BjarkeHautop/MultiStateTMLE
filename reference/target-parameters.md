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
  `~ illness == 0`: only events while it holds count. It may use the
  processes that hazards depend on (their number of events so far).

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
risk(death, illness)
#> Target: risk(death, illness)

# Only deaths before illness count:
risk(death, when = ~ illness == 0)
#> Target: risk(death, when = ~illness == 0)

# For a recurrent process, the risk of a first event rather than the
# expected number of events:
risk(relapse, only_first = TRUE)
#> Target: risk(relapse, only_first = TRUE)

# Processes held in a character vector:
targets <- c("death", "illness")
do.call(risk, as.list(targets))
#> Target: risk(death, illness)

# Use them in tmle_estimate():
fit <- fit_initial(
  illness_death,
  models = list(
    hazard_model(illness ~ biomarker + treatment, limit = 1),
    death ~ biomarker + treatment + illness,
    censoring ~ biomarker
  ),
  tau = 1
)
tmle_estimate(fit, risk(death, illness))
#> TMLE estimate of risk(death, illness)
#> without intervention
#> 
#>     target  estimate         se     lower     upper
#>     <char>     <num>      <num>     <num>     <num>
#> 1:   death 0.3326759 0.03439736 0.2652584 0.4000935
#> 2: illness 0.3687958 0.03488168 0.3004290 0.4371627
#> 
#> death: Initial estimate: 0.3327; no targeting needed, as it already solves the EIC equation
#> illness: Initial estimate: 0.3688; no targeting needed, as it already solves the EIC equation
tmle_estimate(fit, risk(death, when = ~ illness == 0))$estimate
#>    target  estimate         se     lower     upper
#>    <char>     <num>      <num>     <num>     <num>
#> 1:  death 0.1713816 0.02734643 0.1177836 0.2249796
```
