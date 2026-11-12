test_that("fit_initial() fits the illness-death example", {
  fit <- fit_initial(
    illness_death,
    models = illness.death.models,
    tau = 1
  )

  expect_s3_class(fit, "initial_fit")
  expect_equal(fit$n, length(unique(illness_death$id)))
  expect_equal(fit$tau, 1)
  expect_equal(
    fit$processes,
    data.table::data.table(
      name = c("z", "death", "censoring"),
      type = c("one_jump", "terminal", "terminal"),
      model = c("z ~ L0", "death ~ L0 + z", "censoring ~ L0"),
      at_risk = c(TRUE, FALSE, FALSE)
    )
  )
  expect_named(fit$models, c("z", "death", "censoring"))
  expect_s3_class(fit$models$death, "hazard_model")

  expect_equal(fit$process_names, c("z", "death", "censoring"))
  expect_equal(
    fit$process_types,
    list(z = "one_jump", death = "terminal", censoring = "terminal")
  )
  expect_equal(fit$cens_process_id, 3)
  expect_named(fit$at_risks, "z")
  # z is a state variable, with states z = 0 and z = 1.
  expect_equal(fit$depend_matrix$z, c(0, 1))
  expect_true(all(fit$tmp_long$time <= 1))
  expect_true(all(
    c("P.z.1", "P.z.2", "P.death.1", "P.death.2") %in%
      names(fit$tmp_long)
  ))
})

test_that("fit_initial() matches the deprecated fit_types interface", {
  withr::local_options(lifecycle_verbosity = "warning")
  new_fit <- fixture.fit(40)

  expect_warning(
    old_fit <- fit_initial(
      with.delta(fixture.data(40)),
      fit_types = illness.death.fit.types,
      tau = 1
    ),
    class = "lifecycle_warning_deprecated"
  )
  # A fit_types list in place of models is also recognised.
  expect_warning(
    fit_initial(
      with.delta(fixture.data(40)),
      illness.death.fit.types,
      tau = 1
    ),
    class = "lifecycle_warning_deprecated"
  )

  expect_equal(
    tmle_alpha(new_fit, target = "death", alpha = 0.5, z_name = "z")$estimate,
    tmle_alpha(old_fit, target = "death", alpha = 0.5, z_name = "z")$estimate
  )
})

test_that("fit_initial() takes simevent output as is", {
  skip_if_not_installed("simevent")
  sim <- simulate.illness.death(60, seed = 1)

  # simevent's z column, the count after each event, is recomputed.
  expect_true("z" %in% names(sim))
  fit <- fit_initial(sim, models = illness.death.models, tau = 1)
  expect_equal(fit$n, 60)
  expect_equal(fit$depend_matrix$z, c(0, 1))
})

test_that("fit_initial() accepts plain formulas, named lists and any order", {
  fit <- fixture.fit(40)

  reordered <- fit_initial(
    fixture.data(40),
    models = list(
      censoring = censoring ~ L0,
      death = death ~ L0 + z,
      z = hazard_model(z ~ L0, at_risk = ~ z == 0)
    ),
    tau = 1
  )
  expect_equal(
    tmle_alpha(reordered, target = "death", alpha = 0.5, z_name = "z")$estimate,
    tmle_alpha(fit, target = "death", alpha = 0.5, z_name = "z")$estimate
  )
})

test_that("fit_initial() accepts transformed covariates in models", {
  transformed <- illness.death.models
  transformed[[2]] <- death ~ log(exp(L0)) + z

  expect_equal(
    fit_initial(fixture.data(40), models = transformed, tau = 1)$tmp_long,
    fixture.fit(40)$tmp_long
  )
})

test_that("fit_initial() uses a censoring process of another name", {
  data <- fixture.data(40)
  levels(data$event)[levels(data$event) == "censoring"] <- "dropout"

  expect_error(
    fit_initial(
      data,
      models = list(hazard_model(z ~ L0, limit = 1), death ~ L0 + z),
      tau = 1
    ),
    "censoring"
  )
  fit <- fit_initial(
    data,
    models = list(
      hazard_model(z ~ L0, limit = 1),
      death ~ L0 + z,
      dropout ~ L0
    ),
    tau = 1,
    censoring = "dropout"
  )
  expect_equal(fit$process_names[fit$cens_process_id], "dropout")
})

test_that("hazard_model() validates its arguments", {
  expect_error(hazard_model(~L0), "two-sided")
  expect_error(hazard_model(log(z) ~ L0), "two-sided")
  expect_error(hazard_model(z ~ L0, limit = 0), "limit")
  expect_error(hazard_model(z ~ L0, at_risk = z ~ z == 0), "one-sided")
  expect_error(hazard_model(z ~ L0, fit = "glm"), "cox")
  expect_output(
    print(hazard_model(z ~ L0, limit = 1, at_risk = ~ death == 0)),
    "limit: 1"
  )
})

test_that("hazard_model() conditions become at-risk functions", {
  k <- 5
  at_risk <- .hazard_model_to_fit_type(
    hazard_model(z ~ L0, limit = 1, at_risk = ~ death == 0 & z < k),
    delta = 1,
    state_vars = c("z", "death")
  )$at_risk

  # State variables are read from dt, other variables from the formula's
  # environment.
  expect_equal(
    all.names(body(at_risk)),
    c("&", "<", "[[", "dt", "&", "==", "[[", "dt", "<", "[[", "dt", "k")
  )
  dt <- data.table::data.table(z = c(0, 1, 0), death = c(0, 0, 1))
  expect_equal(at_risk(dt), c(TRUE, FALSE, FALSE))
  expect_null(
    .hazard_model_to_fit_type(hazard_model(z ~ L0), 1, "z")$at_risk
  )
})

test_that("fit_initial() validates its arguments", {
  m <- illness.death.models

  expect_error(
    fit_initial(illness_death[, c("id", "time")], models = m, tau = 1),
    "event"
  )
  expect_error(fit_initial(illness_death, models = m, tau = -1), "tau")
  expect_error(fit_initial(illness_death, models = m, tau = 0), "tau")
  expect_error(
    fit_initial(illness_death, models = m[1:2], tau = 1),
    "censoring"
  )
  expect_error(
    fit_initial(illness_death, models = m[c(1, 3)], tau = 1),
    "without a model.*death"
  )
  expect_error(
    fit_initial(illness_death, models = c(m, list(death ~ L0)), tau = 1),
    "one model"
  )
  expect_error(
    fit_initial(illness_death, models = c(m, list("x")), tau = 1),
    "formula or a hazard_model"
  )
  expect_error(
    fit_initial(illness_death, models = list(other = death ~ L0), tau = 1),
    "named \"other\""
  )
  expect_error(
    fit_initial(
      illness_death,
      models = list(
        hazard_model(z ~ L0, at_risk = ~ L0 > 0),
        death ~ L0 + z,
        censoring ~ L0
      ),
      tau = 1
    ),
    "L0"
  )
  data <- illness_death
  data$event[1] <- NA
  expect_error(fit_initial(data, models = m, tau = 1), "missing")
  expect_error(
    fit_initial(illness_death, models = m, tau = 1, arms = 1),
    "treatment_model"
  )
  expect_error(
    fit_initial(
      illness_death,
      models = m,
      tau = 1,
      derived_vars = list(function(dt) dt[["z"]])
    ),
    "derived_vars"
  )
  expect_error(
    fit_initial(
      illness_death,
      models = m,
      tau = 1,
      hal_control = list(foo = 1)
    ),
    "hal_control"
  )
})

test_that("fit_initial() rejects invalid inputs", {
  m <- illness.death.models

  expect_error(fit_initial("x", models = m, tau = 1), "data")
  expect_error(
    fit_initial(fixture.data(), models = m, tau = 1, max_count = -1),
    "max_count"
  )
  expect_error(
    fit_initial(fixture.data(), models = m, tau = 1, verbose = "x"),
    "verbose"
  )
  expect_error(
    fit_initial(
      fixture.data.treated(),
      models = m,
      tau = 1,
      treatment_model = list(model = A0 ~ L0, fit = "glm")
    ),
    "treatment_model"
  )
  expect_error(
    fit_initial(
      fixture.data.treated(),
      models = m,
      tau = 1,
      treatment_model = ~L0,
      arms = 1
    ),
    "two-sided"
  )
})

test_that("fit_initial() converts a non-numeric time column with a warning", {
  data <- fixture.data(20)
  data$time <- as.character(data$time)

  expect_warning(
    fit <- fit_initial(data, models = illness.death.models, tau = 1),
    "not numeric"
  )
  expect_true(is.numeric(fit$tmp_long$time))
})

test_that("fit_initial() fits a treatment model and propensities under 'a'", {
  fit <- fit_initial(
    fixture.data.treated(),
    models = illness.death.models,
    tau = 1,
    treatment_model = A0 ~ L0
  )

  # arms defaults to both arms when a treatment model is given.
  expect_equal(fit$arms, 0:1)
  expect_true(data.table::is.data.table(fit$tmp_long))
  expect_true(all(fit$tmp_long$pi.A0.1 >= 0 & fit$tmp_long$pi.A0.1 <= 1))
  expect_true(all(
    c("clever.weight.a0", "clever.weight.a1") %in%
      names(fit$tmp_long)
  ))
})

test_that("fit_initial() reports fit diagnostics as messages when verbose = TRUE", {
  msgs <- capture_messages(fixture.fit(20, verbose = TRUE))
  expect_true(any(grepl("coxph", msgs)))
  expect_silent(suppressMessages(fixture.fit(20, verbose = TRUE)))
})

test_that("fit_initial()'s HAL setup identifies the process for its own delta", {
  # Regression test for a `name_hal` indexing bug that used to select every
  # process name instead of just the HAL one and error before reaching the
  # HAL fit itself.
  m <- illness.death.models
  m[[1]] <- hazard_model(z ~ L0, limit = 1, fit = "hal")
  assign("fit.hal", function(...) stop("reached fit.hal"), envir = globalenv())
  withr::defer(rm("fit.hal", envir = globalenv()))

  expect_error(
    fit_initial(fixture.data(30), models = m, tau = 1),
    "reached fit.hal"
  )
})

test_that("fit_initial() errors early when HAL is not available", {
  m <- illness.death.models
  m[[1]] <- hazard_model(z ~ L0, limit = 1, fit = "hal")

  expect_error(
    fit_initial(fixture.data(30), models = m, tau = 1),
    "requires a fit.hal() function",
    fixed = TRUE
  )
})

test_that("fit_initial() validates its control settings", {
  m <- illness.death.models

  expect_error(
    fit_initial(illness_death, m, tau = 1, control = list(foo = TRUE)),
    "control"
  )
  expect_error(
    fit_initial(
      illness_death,
      m,
      tau = 1,
      control = list(prune_states = "yes")
    ),
    "control\\$prune_states"
  )
  expect_equal(
    fit_initial(
      fixture.data(30),
      m,
      tau = 1,
      control = list(prune_states = FALSE, use_exponential = FALSE)
    ),
    fixture.fit(30)
  )
})
