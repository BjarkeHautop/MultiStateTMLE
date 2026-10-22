test_that("fit_initial() fits the illness-death example", {
  fit <- fit_initial(
    illness_death,
    fit_types = illness.death.fit.types,
    tau = 1
  )

  expect_s3_class(fit, "initial_fit")
  expect_equal(fit$n, length(unique(illness_death$id)))
  expect_equal(fit$tau, 1)
  expect_equal(
    fit$processes,
    data.table::data.table(
      name = c("z", "death", "censoring"),
      delta = c(2, 1, 0),
      type = c("one_jump", "terminal", "terminal"),
      model = c("~L0", "~L0 + z", "~L0"),
      at_risk = c(TRUE, FALSE, FALSE)
    )
  )

  expect_equal(fit$process_names, c("z", "death", "censoring"))
  expect_equal(fit$process_deltas, c(2, 1, 0))
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

test_that("fit_initial() accepts transformed covariates in models", {
  transformed <- illness.death.fit.types
  transformed$death$model <- ~ log(exp(L0)) + z

  expect_equal(
    fit_initial(illness_death, fit_types = transformed, tau = 1)$tmp_long,
    fit_initial(
      illness_death,
      fit_types = illness.death.fit.types,
      tau = 1
    )$tmp_long
  )
})

test_that("fit_initial() validates its arguments", {
  ft <- illness.death.fit.types

  expect_error(
    fit_initial(illness_death[, c("id", "time")], fit_types = ft, tau = 1),
    "delta"
  )
  expect_error(fit_initial(illness_death, fit_types = ft, tau = -1), "tau")
  expect_error(
    fit_initial(illness_death, fit_types = ft[c("z", "death")], tau = 1),
    "censoring"
  )
  expect_error(
    fit_initial(
      illness_death,
      fit_types = c(ft, list(other = list(model = ~L0))),
      tau = 1
    ),
    "delta"
  )
  expect_error(
    fit_initial(
      illness_death,
      fit_types = c(ft, list(other = list(delta = 1, model = ~L0))),
      tau = 1
    ),
    "own delta code"
  )
  expect_error(
    fit_initial(
      illness_death,
      fit_types = c(ft, list(other = list(delta = 3, model = event ~ L0))),
      tau = 1
    ),
    "one-sided"
  )
  expect_error(
    fit_initial(
      illness_death,
      fit_types = c(ft, list(other = list(delta = 3, model = ~L0, foo = 1))),
      tau = 1
    ),
    "foo"
  )
  expect_error(
    fit_initial(
      illness_death,
      fit_types = list(z = list(delta = 2, model = ~L0, fit = "glm")),
      tau = 1
    ),
    "fit"
  )
  expect_error(
    fit_initial(illness_death, fit_types = ft, tau = 1, arms = 1),
    "treatment_model"
  )
  expect_error(fit_initial(illness_death, fit_types = ft, tau = 0), "tau")
  expect_error(
    fit_initial(
      illness_death,
      fit_types = ft,
      tau = 1,
      derived_vars = list(function(dt) dt[["z"]])
    ),
    "derived_vars"
  )
  expect_error(
    fit_initial(
      illness_death,
      fit_types = ft,
      tau = 1,
      hal_control = list(foo = 1)
    ),
    "hal_control"
  )
})

test_that("fit_initial() rejects invalid inputs", {
  ft <- illness.death.fit.types

  expect_error(fit_initial("x", fit_types = ft, tau = 1), "data")
  expect_error(
    fit_initial(fixture.data(), fit_types = ft, tau = 1, max_count = -1),
    "max_count"
  )
  expect_error(
    fit_initial(fixture.data(), fit_types = ft, tau = 1, verbose = "x"),
    "verbose"
  )
  expect_error(
    fit_initial(
      fixture.data.treated(),
      fit_types = ft,
      tau = 1,
      treatment_model = list(model = A0 ~ L0, fit = "glm")
    ),
    "treatment_model"
  )
  expect_error(
    fit_initial(
      fixture.data.treated(),
      fit_types = ft,
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
    fit <- fit_initial(data, fit_types = illness.death.fit.types, tau = 1),
    "not numeric"
  )
  expect_true(is.numeric(fit$tmp_long$time))
})

test_that("fit_initial() fits a treatment model and propensities under 'a'", {
  fit <- fit_initial(
    fixture.data.treated(),
    fit_types = illness.death.fit.types,
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
  # Regression test for a `name.hal` indexing bug that used to select every
  # process name instead of just the HAL one and error before reaching the
  # HAL fit itself.
  ft <- illness.death.fit.types
  ft$z$fit <- "hal"
  assign("fit.hal", function(...) stop("reached fit.hal"), envir = globalenv())
  withr::defer(rm("fit.hal", envir = globalenv()))

  expect_error(
    fit_initial(fixture.data(30), fit_types = ft, tau = 1),
    "reached fit.hal"
  )
})

test_that("fit_initial() errors early when HAL is not available", {
  ft <- illness.death.fit.types
  ft$z$fit <- "hal"

  expect_error(
    fit_initial(fixture.data(30), fit_types = ft, tau = 1),
    "requires a fit.hal() function",
    fixed = TRUE
  )
})

test_that("fit_initial() validates its control settings", {
  ft <- illness.death.fit.types

  expect_error(
    fit_initial(illness_death, ft, tau = 1, control = list(foo = TRUE)),
    "control"
  )
  expect_error(
    fit_initial(
      illness_death,
      ft,
      tau = 1,
      control = list(prune_states = "yes")
    ),
    "control\\$prune_states"
  )
  expect_equal(
    fit_initial(
      fixture.data(30),
      ft,
      tau = 1,
      control = list(prune_states = FALSE, use_exponential = FALSE)
    ),
    fixture.fit(30)
  )
})
