test_that("as_sim_model() takes types and limits from the fit", {
  skip_if_not_installed("simevent")
  fit <- fixture.fit(80)
  model <- as_sim_model(fit, fixture.data(80))

  expect_named(model$processes, c("illness", "death", "censoring"))
  expect_equal(model$processes$illness$type, "transient")
  expect_equal(model$processes$illness$limit, 1)
  expect_equal(model$processes$death$type, "terminal")
  expect_equal(model$processes$censoring$type, "censoring")

  simulated <- simevent::simulate_events(model, n = 50, seed = 1)
  expect_true(all(table(simulated$id[simulated$event == "illness"]) == 1))
})

test_that("as_sim_model() gives a one-jump process without a limit limit 1", {
  skip_if_not_installed("simevent")
  models <- illness.death.models
  models[[1]] <- illness ~ biomarker
  fit <- fit_initial(fixture.data(80), models = models, tau = 1)

  expect_equal(fit$processes$type[fit$processes$name == "illness"], "one_jump")
  model <- as_sim_model(fit, fixture.data(80))
  expect_equal(model$processes$illness$limit, 1)
})

test_that("as_sim_model() simulates the effects of a recurrent process", {
  skip_if_not_installed("simevent", minimum_version = "0.1.1.9001")
  fit <- recurrent.fit(n = 80, max_count = 2)
  model <- as_sim_model(fit, simulate.recurrent(80, 1))

  expect_equal(model$processes$relapse$limit, Inf)
  froms <- vapply(
    Filter(function(effect) effect$to == "death", model$effects),
    `[[`,
    character(1),
    "from"
  )
  expect_setequal(froms, c("biomarker", "I(relapse >= 1)", "I(relapse >= 2)"))
  expect_s3_class(
    simevent::simulate_events(model, n = 50, seed = 1),
    "data.table"
  )
})

test_that("as_sim_model() rejects what simevent can't simulate", {
  skip_if_not_installed("simevent")
  data <- fixture.data(40)

  at_risk <- illness.death.models
  at_risk[[2]] <- hazard_model(
    death ~ biomarker + illness,
    at_risk = ~ illness == 1
  )
  expect_error(
    as_sim_model(fit_initial(data, models = at_risk, tau = 1), data),
    "`at_risk` condition.*death"
  )

  timed <- fit_initial(
    data,
    models = list(
      hazard_model(illness ~ biomarker, limit = 1),
      death ~ biomarker + illness + T.illness,
      censoring ~ biomarker
    ),
    tau = 1,
    depend_time = list(illness = ~ 1 * (time > 0.5))
  )
  expect_error(as_sim_model(timed, data), "depend_time.*T.illness")
})
