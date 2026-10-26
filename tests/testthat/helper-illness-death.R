# Shared helpers for simulation-based recovery tests: simulate from a
# known illness-death DGP (L0 -> z -> death, z one-jump, death
# terminal, both Weibull(nu=1) i.e. exponential hazards) and compute
# the true P(death by tau) in closed form.
#
# Truth is:
# P(death by tau) = P(direct, from z=0) + P(via a
# jump to z=1 at some u<=tau, then death by tau), integrated over
# L0 ~ N(0,1) numerically.
analytic.truth <- function(
  alpha,
  h_z0 = 0.2,
  h_out0 = 0.3,
  h_out1 = 0.3 * exp(0.7),
  tau = 1
) {
  Q <- function(L0) {
    h_z <- alpha * h_z0 * exp(0.5 * L0)
    p_direct <- (h_out0 / (h_z + h_out0)) * (1 - exp(-(h_z + h_out0) * tau))
    f <- function(u) {
      h_z * exp(-(h_z + h_out0) * u) * (1 - exp(-h_out1 * (tau - u)))
    }
    p_via <- integrate(f, 0, tau)$value
    p_direct + p_via
  }
  integrand <- function(L0) sapply(L0, Q) * dnorm(L0)
  integrate(integrand, -Inf, Inf)$value
}

# Simulate n individuals from the DGP above with simevent, in the
# id/time/delta format of fit_initial().
simulate.illness.death <- function(n, seed) {
  withr::local_seed(seed)
  model <- simevent::sim_model(
    L0 = simevent::sim_covariate(function(N) rnorm(N)),
    z = simevent::sim_process("transient", limit = 1, eta = 0.2, nu = 1),
    death = simevent::sim_process("terminal", eta = 0.3, nu = 1),
    censoring = simevent::sim_process("censoring", eta = 0.1, nu = 1),
    effects = list(
      simevent::sim_effect("L0", "z", 0.5),
      simevent::sim_effect("z", "death", 0.7)
    )
  )
  sim <- simevent::sim_events(model, n = n)
  codes <- c(censoring = 0L, death = 1L, z = 2L)
  data.frame(
    id = sim$id,
    time = sim$time,
    delta = unname(codes[as.character(sim$event)]),
    L0 = sim$L0
  )
}

illness.death.fit.types <- list(
  z = list(
    delta = 2,
    model = ~L0,
    at_risk = function(dt) dt[["z"]] == 0
  ),
  death = list(delta = 1, model = ~ L0 + z),
  censoring = list(delta = 0, model = ~L0)
)

simulate.recovery.fit <- function(alpha, n, seed) {
  initial_fit <- fit_initial(
    simulate.illness.death(n, seed),
    fit_types = illness.death.fit.types,
    tau = 1
  )
  tmle_alpha(initial_fit, target = "death", alpha = alpha, z_name = "z")
}

# The first n individuals of illness_death, for fast unit tests.
fixture.data <- function(n = 40) {
  illness_death[illness_death$id <= n, ]
}

fixture.fit <- function(n = 40, ...) {
  fit_initial(
    fixture.data(n),
    fit_types = illness.death.fit.types,
    tau = 1,
    ...
  )
}

# fixture.data() with a binary baseline treatment A0.
fixture.data.treated <- function(n = 40) {
  data <- fixture.data(n)
  a0 <- withr::with_seed(1, rbinom(n, 1, 0.5))
  data$A0 <- a0[data$id]
  data
}

# Minimal tmle_alpha() and calibrate_alpha() results with the given
# estimates and influence curves, for unit tests of calibrated_contrasts().
# For calibrate_alpha(), `estimate` and `eic` have one entry per target.
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
      estimand = estimand,
      a = a,
      eic = stats::setNames(list(eic), target)
    ),
    class = "tmle_alpha"
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
      a = a,
      estimand = estimand,
      eic = stats::setNames(eic, target)
    ),
    class = "calibrate_alpha"
  )
}

# What solve_alpha() and estimate_derivative() read from a tmle_alpha() fit.
fake.estimate <- function(estimate, se = NA, eic = 0) {
  list(estimate = data.frame(estimate = estimate, se = se), eic = list(eic))
}
