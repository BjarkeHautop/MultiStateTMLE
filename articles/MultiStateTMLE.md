# Getting Started with MultiStateTMLE

``` r

library(MultiStateTMLE)
```

We will in this vignette use the following example: how many deaths
would be prevented by a drug that makes an illness less likely?

## The Data

`illness_death` is a simulated dataset of a trial of 200 patients who
may fall ill and may die, with follow-up subject to censoring. The
`illness` can occur at most once and raises the hazard of death.
Patients were randomised to `treatment` (`1`) or placebo (`0`), and
patients with a higher baseline `biomarker` are more likely to fall ill.

``` r

head(illness_death)
#>   id      time     event  biomarker treatment
#> 1  1 2.8926838 censoring  0.2724785         1
#> 2  2 0.1282264   illness  0.3572619         0
#> 3  2 1.2558891     death  0.3572619         0
#> 4  3 0.6567113     death -0.8616620         0
#> 5  4 0.1752960   illness  0.8083350         1
#> 6  4 0.2409278     death  0.8083350         1
```

Each row is an event: the patient `id`, the `time` of the event, which
`event` it was, and the baseline covariates. A patient who dies without
falling ill has one row, like patient 3, while one who falls ill and
then dies has two, like patient 2:

``` r

illness_death[illness_death$id %in% c(2, 3), ]
#>   id      time   event  biomarker treatment
#> 2  2 0.1282264 illness  0.3572619         0
#> 3  2 1.2558891   death  0.3572619         0
#> 4  3 0.6567113   death -0.8616620         0
```

Simulated data

This is the format
[`simevent::simulate_events()`](https://rdrr.io/pkg/simevent/man/simulate_events.html)
returns, so data simulated with
[simevent](https://github.com/BjarkeHautop/simevent) can be used
directly. `illness_death` itself was simulated with simevent, which lets
us check the estimates against the truth below.

## Fitting the Hazards

The first step is to model the hazard of each event.
[`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md)
takes one formula per event, with the event on the left and the
covariates of its hazard on the right. Using another event’s name as a
covariate gives its number of events so far, so
`death ~ biomarker + treatment + illness` lets the illness change the
death hazard:

``` r

models <- list(
  hazard_model(illness ~ biomarker + treatment, limit = 1),
  death ~ biomarker + treatment + illness,
  censoring ~ biomarker
)
fit <- fit_initial(illness_death, models = models, tau = 1)
fit
#> Initial fit of 200 individuals up to tau = 1
#> 
#>         name     type                                   model at_risk
#>       <char>   <char>                                  <char>  <lgcl>
#> 1:   illness one_jump         illness ~ biomarker + treatment    TRUE
#> 2:     death terminal death ~ biomarker + treatment + illness   FALSE
#> 3: censoring terminal                   censoring ~ biomarker   FALSE
```

By default each hazard is a Cox model, so the death hazard is

\begin{aligned} \lambda\_{\text{death}}(t) = \lambda_0(t) \exp\big( &
\beta_1 \\ \text{biomarker} + \beta_2 \\ \text{treatment} \\ & + \beta_3
\\ \text{illness}(t) \big), \end{aligned}

where \text{illness}(t) is 1 once the patient has fallen ill. `tau = 1`
is the time horizon: everything below concerns what happens by time 1.

Hazard models

Every event in the data needs a model, censoring included. A plain
formula is short for `hazard_model(formula)`, which also takes `limit`,
the most events a patient can have (here one illness), and `at_risk`, a
condition such as `~ illness == 0` for when patients are at risk.

## Halving the Illness

Without any intervention,
[`tmle_estimate()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_estimate.md)
estimates the risk of death by `tau`, given by
[`risk()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/target-parameters.md):

``` r

observed <- tmle_estimate(fit, risk(death))
observed
#> TMLE estimate of risk(death)
#> without intervention
#> 
#>    target  estimate         se     lower     upper
#>    <char>     <num>      <num>     <num>     <num>
#> 1:  death 0.3326759 0.03439736 0.2652584 0.4000935
#> 
#> Initial estimate: 0.3327; no targeting needed, as it already solves the EIC equation
```

Now suppose the drug halves the hazard of the illness.
[`scale_intensity()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/interventions.md)
multiplies a hazard by a factor \alpha, and the estimate is the risk of
death had the illness hazard been \alpha \lambda\_{\text{illness}}(t):

``` r

halved <- tmle_estimate(fit, risk(death), scale_intensity(illness = 0.5))
halved
#> TMLE estimate of risk(death)
#> under scale_intensity(illness = 0.5)
#> 
#>    target  estimate         se     lower     upper
#>    <char>     <num>      <num>     <num>     <num>
#> 1:  death 0.2906665 0.03341164 0.2251809 0.3561521
#> 
#> Initial estimate: 0.2907; no targeting needed, as it already solves the EIC equation
```

Halving the illness lowers the risk of death from 0.333 to 0.291. The
other hazards, death given illness included, are left as they are.

[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`confint()`](https://rdrr.io/r/stats/confint.html) extract the estimate
and its confidence interval:

``` r

coef(halved)
#>     death 
#> 0.2906665
confint(halved, level = 0.9)
#>             5 %      95 %
#> death 0.2357092 0.3456237
```

Reading the output

The standard error comes from the efficient influence curve, so no
bootstrap is needed. The last line gives the estimate from the initial
fit alone, which TMLE then *targets*: it updates the fitted hazards
until the estimate solves the efficient influence curve equation. This
removes the first-order bias of the initial fit, so the confidence
interval is valid even with flexible hazard models.

Here no targeting was needed: the Cox fits already solve the equation to
within the stopping criterion, \text{se} / \log n, so the estimate is
the initial estimate. With 200 patients this is common, since the
remaining bias is small next to the standard error. With larger samples,
or more flexible initial fits, expect one or more targeting steps.

Other interventions

A factor of `0` removes the illness, and one above `1` increases it. The
factor can also be a formula in `time`, the baseline covariates and the
events so far, for example a drug that only works early, or only for
patients with `biomarker > 0`:

``` r

coef(tmle_estimate(
  fit,
  risk(death),
  scale_intensity(illness = ~ ifelse(time < 0.5, 0.5, 1))
))
#>    death 
#> 0.302683
coef(tmle_estimate(
  fit,
  risk(death),
  scale_intensity(illness = ~ ifelse(biomarker > 0, 0.5, 1))
))
#>     death 
#> 0.3053054
```

Several events can be scaled at once, as in
`scale_intensity(illness = 0.5, death = 0.8)`.

## From Hazards to Risks

“Halve the hazard” is hard to interpret: what does it mean for how many
patients fall ill?
[`intervention_curve()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/intervention_curve.md)
answers this for a grid of factors \alpha, estimating both the risk of
the illness, the *exposure*, and the risk of death:

``` r

curve <- intervention_curve(
  fit,
  risk(death),
  exposure = risk(illness),
  alpha = c(0, 0.25, 0.5, 1, 1.5, 2)
)
curve
#> Intervention curve of risk(death) against risk(illness)
#> over 6 scalings of the intensity of illness
#> 
#>    alpha target exposure_estimate exposure_se  estimate         se     lower
#>    <num> <char>             <num>       <num>     <num>      <num>     <num>
#> 1:  0.00  death         0.0000000  0.00000000 0.2336750 0.03966975 0.1559237
#> 2:  0.25  death         0.1145766  0.01453457 0.2644432 0.03553681 0.1947923
#> 3:  0.50  death         0.2169103  0.02428345 0.2906665 0.03374676 0.2245241
#> 4:  1.00  death         0.3687958  0.03488168 0.3326759 0.03427321 0.2655017
#> 5:  1.50  death         0.4790098  0.03999024 0.3645547 0.03712284 0.2917953
#> 6:  2.00  death         0.5616300  0.04251611 0.3893948 0.04055059 0.3099171
#>        upper
#>        <num>
#> 1: 0.3114263
#> 2: 0.3340941
#> 3: 0.3568089
#> 4: 0.3998502
#> 5: 0.4373141
#> 6: 0.4688724
```

Plotting one against the other shows how the risk of death responds to
the risk of illness:

``` r

ce <- curve$estimate

plot(
  ce$exposure_estimate,
  ce$estimate,
  type = "b",
  pch = 19,
  las = 1,
  bty = "l",
  ylim = range(ce$lower, ce$upper),
  xlab = "Risk of illness by time 1",
  ylab = "Risk of death by time 1"
)
arrows(
  ce$exposure_estimate,
  ce$lower,
  ce$exposure_estimate,
  ce$upper,
  angle = 90,
  code = 3,
  length = 0.04
)
```

![Estimated risk of death plotted against the estimated risk of illness
across interventions on the illness hazard, with 95% confidence
intervals. The risk of death increases with the risk of
illness.](MultiStateTMLE_files/figure-html/intervention-plot-1.png)

Removing the illness altogether (the leftmost point) lowers the risk of
death by about 10 percentage points.

Checking against the truth

`illness_death` was simulated from the model below, so simevent can tell
us the true risks under each intervention, by simulating many patients
without censoring and with the illness hazard scaled by `intervene`:

``` r

library(simevent)

true_model <- sim_model(
  biomarker = sim_covariate(function(N) rnorm(N)),
  treatment = sim_covariate(function(N) rbinom(N, 1, 0.5)),
  illness = sim_process("transient", limit = 1, eta = 0.7, nu = 1),
  death = sim_process("terminal", eta = 0.3, nu = 1),
  censoring = sim_process("censoring", eta = 0.1, nu = 1),
  effects = list(
    sim_effect(illness ~ biomarker + treatment, coef = c(0.5, -1)),
    sim_effect(death ~ illness + treatment, coef = c(1.2, -0.3))
  )
)

true_risks <- function(intervene) {
  sim <- simulate_events(
    true_model,
    n = 1e5,
    cens = 0,
    intervene = intervene,
    seed = 1
  )
  risks <- event_risk(sim, true_model, c("illness", "death"), tau = 1)
  setNames(risks$risk, risks$process)
}

check <- ce[, c("alpha", "estimate", "lower", "upper")]
check$truth <- sapply(
  check$alpha,
  function(alpha) true_risks(list(illness = alpha))["death"]
)
check
#>    alpha  estimate     lower     upper   truth
#>    <num>     <num>     <num>     <num>   <num>
#> 1:  0.00 0.2336750 0.1559237 0.3114263 0.22972
#> 2:  0.25 0.2644432 0.1947923 0.3340941 0.25350
#> 3:  0.50 0.2906665 0.2245241 0.3568089 0.27619
#> 4:  1.00 0.3326759 0.2655017 0.3998502 0.31160
#> 5:  1.50 0.3645547 0.2917953 0.4373141 0.34098
#> 6:  2.00 0.3893948 0.3099171 0.4688724 0.36302
```

Every interval covers the truth, and the estimates follow the true curve
closely, removing the illness included.

## Calibrating the Drug

Often the question is about risks directly: what if the drug prevented
20% of illnesses?
[`calibrate()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/interventions.md)
finds the factor \alpha for which the risk of illness reaches a given
value, and estimates the risk of death under it:

``` r

calibrated <- tmle_estimate(
  fit,
  risk(death),
  calibrate(risk(illness), relative_change(0.8))
)
calibrated
#> TMLE estimate of risk(death)
#> under calibrate(risk(illness), relative_change(0.8))
#> 
#>    quantity  estimate         se     lower     upper
#>      <char>     <num>      <num>     <num>     <num>
#> 1:    alpha 0.7443751 0.01232345 0.7202216 0.7685287
#> 2:    theta 0.2950367 0.02790535 0.2403432 0.3497301
#> 
#> Estimates under the calibrated intervention:
#>    target  estimate         se    lower     upper   se_crude
#>    <char>     <num>      <num>    <num>     <num>      <num>
#> 1:  death 0.3127176 0.03330144 0.247448 0.3779872 0.03338458
```

The first table gives the calibrated \alpha and the risk of illness
\theta it reaches, 80% of the observed 0.369. The second gives the risk
of death under it.

Calibration rules

The target risk of illness can be given in three ways:

- `level(theta)`: the risk is `theta`.
- `absolute_change(delta)`: the risk is `delta` plus its observed value.
- `relative_change(rho)`: the risk is `rho` times its observed value.

Because \alpha is itself estimated, the standard error `se` accounts for
its uncertainty. `se_crude` is the standard error if \alpha were known.

[`years_lost()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/target-parameters.md)
gives the expected years of life lost before `tau` instead of the risk,
and can be used wherever
[`risk()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/target-parameters.md)
can:

``` r

coef(tmle_estimate(
  fit,
  years_lost(death),
  calibrate(risk(illness), relative_change(0.8))
))
#>     death 
#> 0.1596098
```

## Explaining a Treatment Effect

So far `treatment` has only been a covariate. Suppose instead we want
its effect on death, and how much of it is because treated patients fall
ill less often. The calibrated intervention above answers this: give
patients the treatment, but calibrate their risk of illness to what it
would have been on placebo. This splits the total effect into two
contrasts:

- *mediated*: treatment versus the calibrated intervention, the part
  that goes through the illness, and
- *calibrated*: the calibrated intervention versus placebo, the rest.

`treatment_model` tells
[`fit_initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/fit_initial.md)
how treatment was assigned, so that estimates can be made with everyone
treated or everyone on placebo:

``` r

trial_fit <- fit_initial(
  illness_death,
  models = models,
  tau = 1,
  treatment_model = treatment ~ biomarker
)
```

`a` sets everyone’s treatment, and `reference_a = 0` calibrates the risk
of illness to that under placebo:

``` r

treated <- tmle_estimate(trial_fit, risk(death), a = 1)
placebo <- tmle_estimate(trial_fit, risk(death), a = 0)
treated_calibrated <- tmle_estimate(
  trial_fit,
  risk(death),
  calibrate(risk(illness), relative_change(1), reference_a = 0),
  a = 1
)
treated_calibrated
#> TMLE estimate of risk(death)
#> under calibrate(risk(illness), relative_change(1), reference_a = 0), with a = 1
#> 
#>    quantity  estimate         se     lower    upper
#>      <char>     <num>      <num>     <num>    <num>
#> 1:    alpha 2.6814303 0.72116544 1.2679720 4.094889
#> 2:    theta 0.4974118 0.05135867 0.3967507 0.598073
#> 
#> Estimates under the calibrated intervention:
#>    target  estimate         se    lower     upper   se_crude
#>    <char>     <num>      <num>    <num>     <num>      <num>
#> 1:  death 0.3112979 0.06357508 0.186693 0.4359028 0.06411261
```

To give treated patients the risk of illness they would have had on
placebo, their illness hazard must be multiplied by about 2.7.
[`calibrated_contrasts()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/calibrated_contrasts.md)
then combines the three fits:

``` r

contrasts <- calibrated_contrasts(treated, placebo, treated_calibrated)
contrasts
#> Decomposition of the effect of a on death into mediated and calibrated contrasts
#> 
#>    target   contrast    estimate         se      lower         upper
#>    <char>     <char>       <num>      <num>      <num>         <num>
#> 1:  death      total -0.18273318 0.06913193 -0.3182293 -0.0472370913
#> 2:  death   mediated -0.06387838 0.03275458 -0.1280762  0.0003194037
#> 3:  death calibrated -0.11885480 0.08063644 -0.2768993  0.0391897229
```

Treatment lowers the risk of death by 18 percentage points, of which 6
go through the illness.

Checking against the truth

With `true_model` from before, simevent gives the true contrasts too:
simulate each arm without censoring, and for the calibrated arm search
for the factor on the illness hazard that matches the risk of illness on
placebo:

``` r

placebo_risks <- true_risks(list(treatment = 0))
treatment_risks <- true_risks(list(treatment = 1))
true_alpha <- uniroot(
  function(alpha) {
    calibrated <- true_risks(list(treatment = 1, illness = alpha))
    calibrated["illness"] - placebo_risks["illness"]
  },
  interval = c(1, 5),
  tol = 0.01
)$root
calibrated_risks <- true_risks(list(treatment = 1, illness = true_alpha))

check <- contrasts$estimate[, c("contrast", "estimate", "lower", "upper")]
check$truth <- c(
  treatment_risks["death"] - placebo_risks["death"],
  treatment_risks["death"] - calibrated_risks["death"],
  calibrated_risks["death"] - placebo_risks["death"]
)
check
#>      contrast    estimate      lower         upper    truth
#>        <char>       <num>      <num>         <num>    <num>
#> 1:      total -0.18273318 -0.3182293 -0.0472370913 -0.13528
#> 2:   mediated -0.06387838 -0.1280762  0.0003194037 -0.05329
#> 3: calibrated -0.11885480 -0.2768993  0.0391897229 -0.08199
```

Every interval covers the truth. With 200 patients the contrasts are
estimated less precisely than the risks, since each compares two arms of
about 100.

## Where Next

- [`risk()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/target-parameters.md)
  and
  [`years_lost()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/target-parameters.md)
  take several events at once, as in `risk(death, illness)`, count only
  events in some states, as in `risk(death, when = ~ illness == 0)`, and
  for recurrent events give the expected number of events, or the risk
  of a first one with `only_first = TRUE`. See
  [`?"target-parameters"`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/target-parameters.md).
- To simulate event history data, and the truth under interventions, see
  [simevent](https://github.com/BjarkeHautop/simevent).
