# Shared helpers for simulation-based recovery tests: simulate from a
# known illness-death DGP (L0 -> z -> outcome1, z one-jump, outcome1
# terminal, both Weibull(nu=1) i.e. exponential hazards) and compute
# the true P(outcome1 by tau) in closed form.
#
# Truth is:
# P(outcome1 by tau) = P(direct, from z=0) + P(via a
# jump to z=1 at some u<=tau, then outcome1 by tau), integrated over
# L0 ~ N(0,1) numerically.
analytic.truth <- function(alpha, h_z0 = 0.2, h_out0 = 0.3, h_out1 = 0.3 * exp(0.7), tau = 1) {
  Q <- function(L0) {
    h_z <- alpha * h_z0 * exp(0.5 * L0)
    p_direct <- (h_out0 / (h_z + h_out0)) * (1 - exp(-(h_z + h_out0) * tau))
    f <- function(u) h_z * exp(-(h_z + h_out0) * u) * (1 - exp(-h_out1 * (tau - u)))
    p_via <- integrate(f, 0, tau)$value
    p_direct + p_via
  }
  integrand <- function(L0) sapply(L0, Q) * dnorm(L0)
  integrate(integrand, -Inf, Inf)$value
}

simulate.recovery.fit <- function(alpha, n, seed) {
  set.seed(seed)
  baseline <- list(L0 = function(N) rnorm(N))
  processes <- list(
    z = list(type = "one.jump", eta = 0.2, nu = 1),
    outcome1 = list(type = "terminal", eta = 0.3, nu = 1),
    censoring = list(type = "censoring", eta = 0.1, nu = 1)
  )
  effects <- list(
    c("L0", "z", 0.5),
    c("L0", "censoring", 0), # prepare.initial doesn't support intercept-only "~1"
    c("z", "outcome1", 0.7)
  )
  tau <- 1

  dt <- sim.generic(
    baseline = baseline,
    processes = processes,
    effects = effects,
    n = n
  )

  initial.fit <- prepare.initial(
    dt,
    tau = tau,
    fit.types = list(
      z = list(
        model = "Surv(tstart, tstop, delta == 2)~L0",
        fit = "cox",
        at.risk = function(dt) (dt[["z"]] == 0)
      ),
      outcome1 = list(model = "Surv(tstart, tstop, delta == 1)~L0+z", fit = "cox"),
      censoring = list(model = "Surv(tstart, tstop, delta == 0)~L0", fit = "cox")
    ),
    verbose = FALSE
  )

  tmle.alpha.fun(
    initial.fit = initial.fit,
    target = "outcome1",
    tau = tau,
    alpha = alpha,
    output.eic = TRUE,
    verbose = FALSE,
    use.cores = 1
  )
}
