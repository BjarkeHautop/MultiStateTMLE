# Getting started with MultiStateTMLE

MultiStateTMLE estimates causal parameters of multistate event processes
under *stochastic intensity interventions*: hypothetical interventions
that scale the intensity of one event process (for example, an illness)
by a factor `alpha`, and ask how the risk of another process (for
example, death) would change. Estimation uses targeted maximum
likelihood estimation (TMLE), which updates an initial fit of the
hazards so that the resulting estimator is efficient and comes with
influence curve based standard errors.

This vignette walks through the main workflow:

1.  Fit the initial hazard models with
    [`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md).
2.  Target a single intervention with
    [`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md).
3.  Trace how the outcome risk depends on the exposure risk with
    [`calibration_curve()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibration_curve.md).
4.  Calibrate the intervention to a target exposure risk with
    [`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md).
5.  Decompose a treatment effect with
    [`calibrated_contrasts()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrated_contrasts.md).

``` r

library(MultiStateTMLE)
```

## The data

The package ships with `illness_death`, simulated event history data for
150 individuals at risk of an illness `z`, which can occur at most once,
and of `death`. The illness raises the hazard of death, and both hazards
depend on a baseline covariate `L0`. Follow-up is subject to censoring.

``` r

head(illness_death)
#>   id      time delta         L0
#> 1  1 1.0836890     1  0.2724785
#> 2  2 0.6539194     1  0.3572619
#> 3  3 0.1212156     1 -0.8616620
#> 4  4 1.1416233     2  0.8083350
#> 5  4 3.6778032     1  0.8083350
#> 6  5 1.0962284     2  0.6127352
```

The data are in long format with one row per event and the columns:

- `id`: individual identifier.
- `time`: event time.
- `delta`: the integer code of the process that jumped, with `0`
  reserved for censoring. Here `1` is death and `2` is illness.
- Baseline covariates, here `L0`.

An individual who dies without falling ill has a single row, while one
who falls ill and later dies has two rows. For example, individual 3
dies, while individual 4 falls ill and then dies:

``` r

illness_death[illness_death$id %in% c(3, 4), ]
#>   id      time delta        L0
#> 3  3 0.1212156     1 -0.861662
#> 4  4 1.1416233     2  0.808335
#> 5  4 3.6778032     1  0.808335
```

## Step 1: Fit the initial hazard models

[`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md)
fits one hazard model per event process, expands the data to the grid of
observed event times up to the horizon `tau`, and builds the discrete
state space that the targeting step works on.

Each process is specified in `fit_types` by its code in the `delta`
column and a one-sided formula of the covariates of its hazard. Other
processes can be used as covariates by name, here the death hazard
depends on whether the individual has fallen ill (`z`). Because the
illness can occur at most once, we say that individuals are only at risk
of `z` while `z == 0`. Exactly one process must be the censoring process
(`delta = 0`).

``` r

fit_types <- list(
  z = list(
    delta = 2,
    model = ~L0,
    at_risk = function(dt) dt[["z"]] == 0
  ),
  death = list(delta = 1, model = ~ L0 + z),
  censoring = list(delta = 0, model = ~L0)
)

fit <- fit_initial(illness_death, fit_types = fit_types, tau = 1)
fit
#> Initial fit of 150 individuals up to tau = 1
#> 
#>         name delta     type   model at_risk
#>       <char> <num>   <char>  <char>  <lgcl>
#> 1:         z     2 one_jump     ~L0    TRUE
#> 2:     death     1 terminal ~L0 + z   FALSE
#> 3: censoring     0 terminal     ~L0   FALSE
```

The `type` of each process is inferred from the data: death and
censoring end follow-up, and the illness occurs at most once.

By default the hazards are fitted with Cox models. Setting `fit = "hal"`
in a process’s entry uses the highly adaptive lasso instead
(experimental).

## Step 2: Target a single intervention

[`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md)
estimates the risk of the `target` process by time `tau` under the
intervention that multiplies the intensity of `z_name` by `alpha`. With
`alpha = 1` there is no intervention, so this is the observed risk of
death by time 1:

``` r

observed <- tmle_alpha(fit, target = "death", alpha = 1)
observed
#> TMLE estimate without intervention
#> 
#>    target  estimate         se     lower     upper
#>    <char>     <num>      <num>     <num>     <num>
#> 1:  death 0.3250469 0.03906316 0.2484845 0.4016093
#> 
#> Initial estimate: 0.325; converged after 1 iterations
```

The printed table gives the targeted estimate with its standard error
and a 95% confidence interval, and it is stored as `observed$estimate`.
Halving the intensity of the illness gives:

``` r

halved <- tmle_alpha(fit, target = "death", alpha = 0.5, z_name = "z")
halved
#> TMLE estimate under alpha = 0.5 on z
#> 
#>    target  estimate         se     lower     upper
#>    <char>     <num>      <num>     <num>     <num>
#> 1:  death 0.3040549 0.03853505 0.2285276 0.3795823
#> 
#> Initial estimate: 0.3124; converged after 2 iterations
```

[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`confint()`](https://rdrr.io/r/stats/confint.html) extract the estimate
and confidence intervals at other levels:

``` r

coef(halved)
#>     death 
#> 0.3040549
confint(halved, level = 0.9)
#>             5 %      95 %
#> death 0.2406704 0.3674395
```

Setting `alpha = 0` removes the illness altogether, while `alpha > 1`
increases its intensity.

## Step 3: Trace the calibration curve

An intensity scaling `alpha` is not always easy to interpret. It is
often more natural to ask what happens to the outcome when the *risk* of
the exposure is changed by a given amount.
[`calibration_curve()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibration_curve.md)
evaluates a grid of interventions and, for each, estimates both the risk
of the exposure `z` and the risk of the target:

``` r

curve <- calibration_curve(
  fit,
  target = "death",
  z_name = "z",
  alpha_grid = c(0, 0.25, 0.5, 1, 1.5, 2)
)
curve$estimate
#>    alpha target exposure_estimate exposure_se  estimate         se     lower
#>    <num> <char>             <num>       <num>     <num>      <num>     <num>
#> 1:  0.00  death        0.00000000 0.000000000 0.2820265 0.03975707 0.2041041
#> 2:  0.25  death        0.04391205 0.009617945 0.2933270 0.03919098 0.2165141
#> 3:  0.50  death        0.08486244 0.017876912 0.3040549 0.03899060 0.2276348
#> 4:  1.00  death        0.15891283 0.031321396 0.3250469 0.03915660 0.2483014
#> 5:  1.50  death        0.23212809 0.041662682 0.3367827 0.04049368 0.2574166
#> 6:  2.00  death        0.29250854 0.049801208 0.3584115 0.04287637 0.2743753
#>        upper
#>        <num>
#> 1: 0.3599489
#> 2: 0.3701399
#> 3: 0.3804751
#> 4: 0.4017924
#> 5: 0.4161489
#> 6: 0.4424476
```

Plotting the target risk against the exposure risk traces out the
calibration curve:

``` r

ce <- curve$estimate

plot(
  ce$exposure_estimate,
  ce$estimate,
  type = "b",
  pch = 19,
  ylim = range(ce$lower, ce$upper),
  xlab = "Risk of illness by time 1",
  ylab = "Risk of death by time 1"
)
arrows(ce$exposure_estimate, ce$lower, ce$exposure_estimate, ce$upper,
  angle = 90, code = 3, length = 0.04
)
```

![Estimated risk of death plotted against the estimated risk of illness
across interventions on the illness intensity, with 95% confidence
bands. The risk of death increases with the risk of
illness.](MultiStateTMLE_files/figure-html/calibration-plot-1.png)

## Step 4: Calibrate the intervention

[`calibrate_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrate_alpha.md)
solves for the `alpha` that achieves a given exposure risk `theta` and
estimates the target risk under that intervention. The calibration
target is given by a `value` on a `scale`:

- `"absolute"`: `theta = value`.
- `"ratio"`: `theta` is `value` times the observed risk of `z_name`.
- `"difference"`: `theta` is `value` plus the observed risk of `z_name`.

For example, the risk of death if the risk of illness were reduced to
80% of its observed value:

``` r

cal <- calibrate_alpha(
  fit,
  target = "death",
  z_name = "z",
  value = 0.8,
  scale = "ratio"
)
cal
#> Calibrated intervention on z
#> 
#>    quantity  estimate         se      lower     upper
#>      <char>     <num>      <num>      <num>     <num>
#> 1:    alpha 0.7927238 0.00995552 0.77321133 0.8122363
#> 2:    theta 0.1271303 0.02505712 0.07801922 0.1762413
#> 
#> Estimates under the calibrated intervention:
#>    target  estimate         se     lower     upper   se_crude
#>    <char>     <num>      <num>     <num>     <num>      <num>
#> 1:  death 0.3199203 0.03864198 0.2441835 0.3956572 0.03866085
```

The first table gives the calibrated intervention `alpha` and the
exposure risk `theta` it achieves, and the second the risk of death
under it. Its standard error `se` accounts for `alpha` being estimated,
whereas `se_crude` treats it as known.

The target can also be the expected years of life lost before `tau`
rather than the risk of death, with `estimand = "years_lost"`. The
exposure is still calibrated on its risk, unless `exposure_estimand`
says otherwise:

``` r

calibrate_alpha(
  fit,
  target = "death",
  z_name = "z",
  value = 0.8,
  scale = "ratio",
  estimand = "years_lost"
)$estimate
#>    target  estimate         se     lower     upper   se_crude
#>    <char>     <num>      <num>     <num>     <num>      <num>
#> 1:  death 0.1644942 0.02295448 0.1195042 0.2094842 0.02295436
```

## Step 5: Decompose a treatment effect

With a baseline treatment `A0`, the total effect of treatment on the
outcome can be split into the part that goes through its effect on the
exposure and the remainder. The calibrated intervention sets the treated
individuals’ exposure risk to what it would have been under placebo, and
[`calibrated_contrasts()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrated_contrasts.md)
then decomposes the total effect into:

- a *mediated* contrast: treatment versus the calibrated intervention,
  and
- a *calibrated* contrast: the calibrated intervention versus placebo.

`illness_death` has no treatment, so for illustration we add a
randomised treatment that has no actual effect. The treatment model is
given in `treatment_model`, and
[`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md)
then prepares estimates under both treatment values, `a = 0` and
`a = 1`.

``` r

set.seed(1)
treated <- illness_death
a0 <- rbinom(max(treated$id), 1, 0.5)
treated$A0 <- a0[treated$id]

fit_treated <- fit_initial(
  treated,
  fit_types = list(
    z = list(
      delta = 2,
      model = ~ L0 + A0,
      at_risk = function(dt) dt[["z"]] == 0
    ),
    death = list(delta = 1, model = ~ L0 + A0 + z),
    censoring = list(delta = 0, model = ~L0)
  ),
  tau = 1,
  treatment_model = A0 ~ L0
)
```

The contrasts are computed from the influence curves that each fit
returns. In the calibrated arm we set `a = 1` and calibrate to a ratio
of 1 to the risk of illness under placebo, `reference_a = 0`, so the
risk of illness under treatment is calibrated to equal that under
placebo:

``` r

treatment_fit <- tmle_alpha(fit_treated, target = "death", a = 1)
placebo_fit <- tmle_alpha(fit_treated, target = "death", a = 0)
calibrated_fit <- calibrate_alpha(
  fit_treated,
  target = "death",
  z_name = "z",
  value = 1,
  scale = "ratio",
  a = 1,
  reference_a = 0
)

calibrated_contrasts(treatment_fit, placebo_fit, calibrated_fit)
#> Decomposition of the effect of a on death into mediated and calibrated contrasts
#> 
#>    target   contrast   estimate         se       lower      upper
#>    <char>     <char>      <num>      <num>       <num>      <num>
#> 1:  death      total -0.1405813 0.07724416 -0.29197704 0.01081449
#> 2:  death   mediated -0.0234265 0.01632300 -0.05541899 0.00856599
#> 3:  death calibrated -0.1171548 0.08141345 -0.27672220 0.04241266
```

The treatment was randomised without any effect, so the confidence
intervals of all three contrasts cover zero, as expected.

## Further reading

- The function reference documents further options, such as recurrent
  targets, only counting first events (`only_first`), and time-varying
  interventions (`alpha_list`).
- To simulate multistate event history data, see the
  [simevent](https://github.com/BjarkeHautop/simevent) package.
