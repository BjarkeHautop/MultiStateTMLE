# Define a state variable

Describes a variable that
[`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md)
tracks over the event history of each individual, for use in the hazard
models. It starts at `start`, and each event of a process named in `...`
updates it by that process's rule. Pass it to
[`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md)
in `state_variables`, named by the variable, which can then be used by
name in the `models` formulas and their `at_risk` conditions, and in the
`when` condition of a target in
[`tmle_estimate()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_estimate.md).

## Usage

``` r
state_variable(start, ...)
```

## Arguments

- start:

  Number. The value of the variable before any event.

- ...:

  The rules of the variable, named by the process whose events they
  apply to: a number, the value after an event of the process, or a
  vectorised function of the value before the event that returns the
  value after it, such as `\(s) pmin(s + 1, 3)`. Events of other
  processes leave the variable unchanged.

## Value

An object of class `state_variable`.

## Details

A state variable is a state of its own, tracked exactly. Use one for a
quantity of several processes, such as whether treatment is stopped
(stops minus restarts): in a model it would keep changing with the
counts, which the state space of
[`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md)
only tracks up to the count from which on the models stop changing (see
its Counts section).

## Examples

``` r
# Whether treatment is stopped: stop events (stop) set it to 1, and
# restart events (restart) set it back to 0:
state_variable(start = 0, stop = 1, restart = 0)
#> State variable starting at 0
#>   after stop: 1
#>   after restart: 0

# The number of stops, counted up to 3:
state_variable(start = 0, stop = \(s) pmin(s + 1, 3))
#> State variable starting at 0
#>   after stop: function (s) pmin(s + 1, 3)

# \donttest{
# Treatment that is stopped and restarted, any number of times, with a
# death hazard that is higher while off treatment:
if (requireNamespace("simevent", quietly = TRUE)) {
  model <- simevent::sim_model(
    biomarker = simevent::sim_covariate(function(N) rnorm(N)),
    stop = simevent::sim_process("transient", eta = 0.8, nu = 1),
    restart = simevent::sim_process("transient", eta = 0.8, nu = 1),
    death = simevent::sim_process("terminal", eta = 0.2, nu = 1),
    censoring = simevent::sim_process("censoring", eta = 0.2, nu = 1),
    effects = list(
      simevent::sim_effect(stop ~ biomarker, coef = 0.5),
      simevent::sim_effect(death ~ (stop - restart), coef = 1),
      # Stops only while on treatment, restarts only while off it:
      simevent::sim_effect(stop ~ (stop != restart), coef = -1e5),
      simevent::sim_effect(restart ~ (restart != stop - 1), coef = -1e5)
    )
  )
  treatment <- simevent::simulate_events(model, n = 200, seed = 1)

  # Some individuals stop several times. Off treatment as I(stop - restart)
  # keeps changing with the counts, which the state space can't track,
  # while the state variable off is exact:
  fit <- fit_initial(
    treatment,
    models = list(
      hazard_model(stop ~ biomarker, at_risk = ~ off == 0),
      hazard_model(restart ~ biomarker, at_risk = ~ off == 1),
      death ~ biomarker + off,
      censoring ~ biomarker
    ),
    tau = 2,
    state_variables = list(
      off = state_variable(start = 0, stop = 1, restart = 0)
    )
  )
  print(fit)

  # A stop leads to a state with off = 1, and a restart to one with
  # off = 0:
  print(fit$depend_matrix)

  tmle_estimate(fit, risk(death))
}
#> Initial fit of 200 individuals up to tau = 2
#> 
#>         name      type                   model at_risk tracked
#>       <char>    <char>                  <char>  <lgcl>   <int>
#> 1:      stop recurrent        stop ~ biomarker    TRUE       0
#> 2:   restart recurrent     restart ~ biomarker    TRUE       0
#> 3:     death  terminal death ~ biomarker + off   FALSE      NA
#> 4: censoring  terminal   censoring ~ biomarker   FALSE      NA
#> 
#> off: State variable starting at 0
#>   after stop: 1
#>   after restart: 0
#>     stop restart   off state s.stop s.restart
#>    <int>   <int> <num> <int>  <int>     <int>
#> 1:     0       0     0     1      2         1
#> 2:     0       0     1     2      2         1
#> TMLE estimate of risk(death)
#> without intervention
#> 
#>    target  estimate         se     lower     upper
#>    <char>     <num>      <num>     <num>     <num>
#> 1:  death 0.4330546 0.03952905 0.3555791 0.5105301
#> 
#> Initial estimate: 0.4331; no targeting needed, as it already solves the EIC equation
# }
```
