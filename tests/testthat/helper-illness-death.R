# Shared helpers for simulation-based recovery tests: simulate from a
# known illness-death DGP (biomarker -> illness -> death, illness
# one-jump, death terminal, both Weibull(nu=1) i.e. exponential hazards)
# and compute
# the true P(death by tau) in closed form.
#
# Truth is:
# P(death by tau) = P(direct, from illness=0) + P(via a
# jump to illness=1 at some u<=tau, then death by tau), integrated over
# biomarker ~ N(0,1) numerically.
analytic.truth <- function(
  alpha,
  h_z0 = 0.2,
  h_out0 = 0.3,
  h_out1 = 0.3 * exp(0.7),
  tau = 1
) {
  Q <- function(biomarker) {
    h_z <- alpha * h_z0 * exp(0.5 * biomarker)
    p_direct <- (h_out0 / (h_z + h_out0)) * (1 - exp(-(h_z + h_out0) * tau))
    f <- function(u) {
      h_z * exp(-(h_z + h_out0) * u) * (1 - exp(-h_out1 * (tau - u)))
    }
    p_via <- integrate(f, 0, tau)$value
    p_direct + p_via
  }
  integrand <- function(biomarker) sapply(biomarker, Q) * dnorm(biomarker)
  integrate(integrand, -Inf, Inf)$value
}

# Simulate n individuals from the DGP above with simevent, whose output
# fit_initial() takes as is.
simulate.illness.death <- function(n, seed) {
  withr::local_seed(seed)
  model <- simevent::sim_model(
    biomarker = simevent::sim_covariate(function(N) rnorm(N)),
    illness = simevent::sim_process("transient", limit = 1, eta = 0.2, nu = 1),
    death = simevent::sim_process("terminal", eta = 0.3, nu = 1),
    censoring = simevent::sim_process("censoring", eta = 0.1, nu = 1),
    effects = list(
      simevent::sim_effect(illness ~ biomarker, coef = 0.5),
      simevent::sim_effect(death ~ illness, coef = 0.7)
    )
  )
  simevent::simulate_events(model, n = n)
}

illness.death.models <- list(
  hazard_model(illness ~ biomarker, limit = 1),
  death ~ biomarker + illness,
  censoring ~ biomarker
)

simulate.recovery.fit <- function(alpha, n, seed) {
  initial_fit <- fit_initial(
    simulate.illness.death(n, seed),
    models = illness.death.models,
    tau = 1
  )
  tmle_estimate(initial_fit, risk(death), scale_intensity(illness = alpha))
}

# The first n individuals of illness_death, for fast unit tests. Includes
# the randomised baseline treatment A0.
fixture.data <- function(n = 40) {
  illness_death[illness_death$id <= n, ]
}

fixture.fit <- function(n = 40, ...) {
  fit_initial(
    fixture.data(n),
    models = illness.death.models,
    tau = 1,
    ...
  )
}

# Minimal tmle_estimate() results with the given estimates and influence
# curves, for unit tests of calibrated_contrasts(). For fake.calibrated.fit(),
# `estimate` and `eic` have one entry per target.
fake.tmle.fit <- function(
  estimate,
  eic,
  target = "death",
  a = 1,
  estimand = "risk"
) {
  structure(
    list(
      estimate = data.table::data.table(target = target, estimate = estimate),
      target = do.call(estimand, as.list(target)),
      a = a,
      eic = stats::setNames(list(eic), target)
    ),
    class = "tmle_estimate"
  )
}

fake.calibrated.fit <- function(
  estimate,
  eic,
  target = "death",
  a = 1,
  estimand = "risk"
) {
  if (!is.list(eic)) {
    eic <- list(eic)
  }
  structure(
    list(
      estimate = data.table::data.table(target = target, estimate = estimate),
      target = do.call(estimand, as.list(target)),
      intervention = calibrate(risk(illness), relative_change(1)),
      a = a,
      eic = stats::setNames(eic, target)
    ),
    class = "tmle_estimate"
  )
}

# What solve_alpha() and estimate_derivative() read from a .tmle_fixed() fit.
fake.estimate <- function(estimate, se = NA, eic = 0) {
  list(estimate = data.frame(estimate = estimate, se = se), eic = list(eic))
}
