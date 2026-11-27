test_that("state_variable() checks its start and rules", {
  variable <- state_variable(start = 0, stop = 1, restart = \(s) s - 1)
  expect_s3_class(variable, "state_variable")
  expect_equal(variable$start, 0)
  expect_named(variable$rules, c("stop", "restart"))

  expect_error(state_variable(start = 0), "at least one process")
  expect_error(state_variable(start = 0, 1), "named by its process")
  expect_error(
    do.call(
      state_variable,
      stats::setNames(list(0, 1, 0), c("start", "a", "a"))
    ),
    "one rule per process"
  )
  expect_error(state_variable(start = 0, a = "x"), "number or a function")
  expect_error(state_variable(start = NA, a = 1))
  expect_output(print(variable), "State variable starting at 0")
  expect_output(print(variable), "after restart: function")
})

test_that("a state variable gives the same estimate as the process it tracks", {
  ill <- list(ill = state_variable(start = 0, illness = 1))
  fit <- fit_initial(
    fixture.data(60),
    models = list(
      hazard_model(illness ~ biomarker, at_risk = ~ ill == 0),
      death ~ biomarker + ill,
      censoring ~ biomarker
    ),
    tau = 1,
    state_variables = ill
  )
  reference <- fixture.fit(60)

  expect_equal(fit$tmp_long$ill, fit$tmp_long$illness)
  expect_equal(unique(fit$depend_matrix$ill), c(0, 1))
  expect_equal(
    unname(coef(fit$cox_fits$death)),
    unname(coef(reference$cox_fits$death))
  )
  expect_equal(
    tmle_estimate(fit, risk(death))$estimate,
    tmle_estimate(reference, risk(death))$estimate
  )
  expect_output(print(fit), "ill: State variable starting at 0")
})

test_that("a state variable can count events with a cap of its own", {
  skip_if_not_installed("simevent")
  data <- simulate.recurrent(80, seed = 1)
  capped <- list(n = state_variable(start = 0, relapse = \(s) pmin(s + 1, 2)))
  fit <- fit_initial(
    data,
    models = list(
      relapse ~ biomarker,
      death ~ biomarker + n,
      censoring ~ biomarker
    ),
    tau = 1,
    state_variables = capped
  )
  expect_equal(fit$tmp_long$n, pmin(fit$tmp_long$relapse, 2))
  expect_equal(
    sort(unique(fit$depend_matrix$n)),
    sort(unique(fit$tmp_long$n))
  )
  expect_equal(sort(unique(fit$depend_matrix$n)), c(0, 1, 2))
  # No model uses the relapse count itself, so it is not tracked beyond 0,
  # and a relapse moves n up by one, to its cap.
  expect_equal(fit$processes$tracked[[1]], 0L)
  dm <- fit$depend_matrix
  expect_equal(dm$relapse, c(0L, 0L, 0L))
  after <- dm[match(dm$s.relapse, dm$state)]
  expect_equal(after$n, pmin(dm$n + 1, 2))

  uncapped <- list(n = state_variable(start = 0, relapse = \(s) s + 1))
  expect_error(
    fit_initial(
      data,
      models = list(
        relapse ~ biomarker,
        death ~ biomarker + n,
        censoring ~ biomarker
      ),
      tau = 1,
      state_variables = uncapped
    ),
    "does not take in the data up to tau"
  )
})

test_that("fit_initial() checks the state variables", {
  fit_with <- function(state_variables, models = illness.death.models) {
    fit_initial(
      fixture.data(40),
      models = models,
      tau = 1,
      state_variables = state_variables
    )
  }
  expect_error(
    fit_with(list(illness = state_variable(0, illness = 1))),
    "name of a process"
  )
  expect_error(
    fit_with(list(biomarker = state_variable(0, illness = 1))),
    "name of a column of data"
  )
  expect_error(
    fit_with(list(state = state_variable(0, illness = 1))),
    "used internally"
  )
  expect_error(
    fit_with(list(ill = state_variable(0, relapse = 1))),
    "relapse, which have no model"
  )
  expect_error(
    fit_with(list(ill = state_variable(0, censoring = 1))),
    "censoring process"
  )
  expect_error(fit_with(list(state_variable(0, illness = 1))))
  expect_error(
    fit_with(
      list(ill = state_variable(0, illness = \(s) "a")),
      models = list(
        hazard_model(illness ~ biomarker, limit = 1),
        death ~ biomarker + ill,
        censoring ~ biomarker
      )
    ),
    "must return a finite number"
  )
})

test_that("as_sim_model() refuses fits with state variables", {
  skip_if_not_installed("simevent")
  fit <- fit_initial(
    fixture.data(40),
    models = list(
      hazard_model(illness ~ biomarker, limit = 1),
      death ~ biomarker + ill,
      censoring ~ biomarker
    ),
    tau = 1,
    state_variables = list(ill = state_variable(0, illness = 1))
  )
  expect_error(as_sim_model(fit, fixture.data(40)), "state variables")
})

# Treatment that is stopped (Nz) and restarted (Nl), alternately: whether it
# is stopped is Nz - Nl, which keeps changing with the counts.
test_that("a stopped-treatment state variable tracks what the counts can't", {
  skip_if_not_installed("simevent")
  big <- 1e5
  model <- simevent::sim_model(
    L0 = simevent::sim_covariate(function(N) runif(N)),
    censoring = simevent::sim_process("censoring", eta = 0.1, nu = 1.1),
    outcome1 = simevent::sim_process("terminal", eta = 0.025, nu = 1.1),
    Nz = simevent::sim_process("transient", eta = 0.3, nu = 1.1),
    Nl = simevent::sim_process("transient", eta = 0.12, nu = 1.1),
    Nw = simevent::sim_process("transient", eta = 0.1, nu = 1.1, limit = 1),
    cr1 = simevent::sim_process("terminal", eta = 0.025, nu = 1.1),
    effects = list(
      simevent::sim_effect(Nz ~ L0, coef = 1),
      simevent::sim_effect(Nl ~ L0, coef = 1),
      simevent::sim_effect(Nw ~ L0, coef = 1),
      simevent::sim_effect(Nz ~ Nw, coef = 1.5),
      simevent::sim_effect(outcome1 ~ Nw, coef = 2),
      simevent::sim_effect(outcome1 ~ (Nz - Nl), coef = 1),
      simevent::sim_effect(Nz ~ (Nz != Nl), coef = -big),
      simevent::sim_effect(Nl ~ (Nl != Nz - 1), coef = -big)
    )
  )
  data <- simevent::simulate_events(model, n = 150, max_cens = 5, seed = 1)
  data$event <- ifelse(
    data$event == "max_cens",
    "censoring",
    as.character(data$event)
  )
  expect_gt(max(table(data$id[data$event == "Nz"])), 2)

  # The effects of the first two stops and restarts.
  counts_z <- ~ L0 + I(Nl >= 1) + I(Nl >= 2) + I(Nz >= 1) + I(Nz >= 2)
  counts_zw <- stats::update(counts_z, ~ . + Nw)
  fit_with <- function(outcome, at_risk_z, at_risk_l, ...) {
    suppressWarnings(fit_initial(
      data,
      models = list(
        hazard_model(stats::update(counts_zw, Nz ~ .), at_risk = at_risk_z),
        hazard_model(stats::update(counts_zw, Nl ~ .), at_risk = at_risk_l),
        hazard_model(stats::update(counts_z, Nw ~ .), limit = 1),
        outcome,
        stats::update(counts_zw, cr1 ~ .),
        stats::update(counts_zw, censoring ~ .)
      ),
      tau = 3,
      ...
    ))
  }
  off <- fit_with(
    outcome1 ~ L0 + Nw + off,
    ~ off == 0,
    ~ off == 1,
    state_variables = list(off = state_variable(start = 0, Nz = 1, Nl = 0))
  )

  expect_equal(off$processes$tracked[1:3], c(2L, 2L, 1L))
  expect_equal(off$tmp_long$off, off$tmp_long$Nz - off$tmp_long$Nl)
  dm <- off$depend_matrix
  expect_equal(sort(unique(dm$off)), c(0, 1))
  # Only the states reachable by alternating stops and restarts are kept,
  # and a stop sets off to 1, a restart to 0, whatever the counts.
  expect_true(all(dm[off == 0]$Nz == dm[off == 0]$Nl | dm[off == 0]$Nz == 2))
  expect_equal(
    dm[match(dm$s.Nz, dm$state)][dm$off == 0]$off,
    rep(1, sum(dm$off == 0))
  )
  expect_equal(
    dm[match(dm$s.Nl, dm$state)][dm$off == 1]$off,
    rep(0, sum(dm$off == 1))
  )
  expect_true(tmle_estimate(off, risk(outcome1))$checks$converged)

  # From the counts, off treatment keeps changing with them.
  expect_error(
    fit_with(outcome1 ~ L0 + Nw + I(Nz - Nl), ~ off == 0, ~ off == 1),
    "define a state_variable"
  )
  expect_error(
    fit_with(outcome1 ~ L0 + Nw, ~ Nz == Nl, ~ Nl == Nz - 1),
    "The at-risk condition of Nz keeps changing"
  )
})
