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

test_that("a state variable is tracked beyond max_count", {
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
    max_count = 1,
    state_variables = capped
  )
  expect_equal(fit$tmp_long$n, pmin(fit$tmp_long$relapse, 2))
  expect_equal(
    sort(unique(fit$depend_matrix$n)),
    sort(unique(fit$tmp_long$n))
  )
  expect_equal(sort(unique(fit$depend_matrix$n)), c(0, 1, 2))
  # The relapse count is capped at max_count, the state variable is not:
  # after a relapse, n goes up by one in every state below the cap.
  dm <- fit$depend_matrix
  after <- dm[match(dm$s.relapse, dm$state)]
  expect_equal(after$n, pmin(dm$n + 1, max(dm$n)))
  expect_equal(after$relapse, pmin(dm$relapse + 1, 1))

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
# is stopped is Nz - Nl, which the capped counts only give up to max_count.
test_that("a stopped-treatment state variable does not depend on max_count", {
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
  max_stops <- max(table(data$id[data$event == "Nz"]))
  expect_gt(max_stops, 2)

  # The effects of the first two stops and restarts, whatever max_count.
  counts_z <- ~ L0 + I(Nl >= 1) + I(Nl >= 2) + I(Nz >= 1) + I(Nz >= 2)
  counts_zw <- stats::update(counts_z, ~ . + Nw)
  fit_off <- function(max_count, prune_states) {
    suppressWarnings(fit_initial(
      data,
      models = list(
        hazard_model(stats::update(counts_zw, Nz ~ .), at_risk = ~ off == 0),
        hazard_model(stats::update(counts_zw, Nl ~ .), at_risk = ~ off == 1),
        hazard_model(stats::update(counts_z, Nw ~ .), limit = 1),
        outcome1 ~ L0 + Nw + off,
        stats::update(counts_zw, cr1 ~ .),
        stats::update(counts_zw, censoring ~ .)
      ),
      tau = 3,
      max_count = max_count,
      control = list(prune_states = prune_states),
      state_variables = list(off = state_variable(start = 0, Nz = 1, Nl = 0))
    ))
  }
  off <- fit_off(2, prune_states = TRUE)

  expect_equal(off$tmp_long$off, off$tmp_long$Nz - off$tmp_long$Nl)
  dm <- off$depend_matrix
  expect_equal(sort(unique(dm$off)), c(0, 1))
  expect_equal(dm[match(dm$s.Nz, dm$state)]$off, rep(1, nrow(dm)))
  expect_equal(dm[match(dm$s.Nl, dm$state)]$off, rep(0, nrow(dm)))

  estimate <- tmle_estimate(off, risk(outcome1))$estimate
  expect_equal(
    tmle_estimate(
      fit_off(max_stops, prune_states = FALSE),
      risk(outcome1)
    )$estimate,
    estimate
  )

  # With I(Nz - Nl) and every stop tracked, the counts are capped at the
  # largest counts in the data up to tau, so a jump from there leaves Nz - Nl
  # unchanged; the estimate is close but not equal.
  every_count <- suppressWarnings(fit_initial(
    data,
    models = list(
      hazard_model(stats::update(counts_zw, Nz ~ .), at_risk = ~ Nz - Nl == 0),
      hazard_model(stats::update(counts_zw, Nl ~ .), at_risk = ~ Nz - Nl == 1),
      hazard_model(stats::update(counts_z, Nw ~ .), limit = 1),
      outcome1 ~ L0 + Nw + I(Nz - Nl),
      stats::update(counts_zw, cr1 ~ .),
      stats::update(counts_zw, censoring ~ .)
    ),
    tau = 3,
    max_count = max_stops,
    control = list(prune_states = TRUE)
  ))
  for (name in names(off$cox_fits)) {
    expect_equal(
      off$cox_fits[[name]]$loglik,
      every_count$cox_fits[[name]]$loglik
    )
  }
  expect_equal(
    tmle_estimate(every_count, risk(outcome1))$estimate,
    estimate,
    tolerance = 0.01
  )
})
