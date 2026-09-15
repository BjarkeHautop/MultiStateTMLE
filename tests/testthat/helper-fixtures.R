# Shared small illness-death fixtures for unit/validation tests that don't
# need the large-n recovery setup in helper-recovery.R.

make.fixture.dt <- function(n = 40, seed = 1405) {
  set.seed(seed)
  baseline <- list(L0 = function(N) rnorm(N))
  processes <- list(
    z = list(type = "one.jump", eta = 0.2, nu = 1),
    outcome1 = list(type = "terminal", eta = 0.3, nu = 1),
    censoring = list(type = "censoring", eta = 0.1, nu = 1)
  )
  effects <- list(
    c("L0", "z", 0.5),
    c("z", "outcome1", 0.7)
  )
  sim.generic(
    baseline = baseline,
    processes = processes,
    effects = effects,
    n = n
  )
}

make.fixture.fit.types <- function() {
  list(
    z = list(
      model = "Surv(tstart, tstop, delta == 2)~L0",
      fit = "cox",
      at.risk = function(dt) (dt[["z"]] == 0)
    ),
    outcome1 = list(model = "Surv(tstart, tstop, delta == 1)~L0+z", fit = "cox"),
    censoring = list(model = "Surv(tstart, tstop, delta == 0)~L0", fit = "cox")
  )
}

make.fixture.initial.fit <- function(n = 40, seed = 1405, tau = 1, ...) {
  dt <- make.fixture.dt(n = n, seed = seed)
  prepare.initial(
    dt,
    tau = tau,
    fit.types = make.fixture.fit.types(),
    verbose = FALSE,
    ...
  )
}
