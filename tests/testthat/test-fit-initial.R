test_that("fit_initial() fits the illness-death example", {
  fit <- fit_initial(
    illness_death,
    fit_types = illness.death.fit.types,
    tau = 1
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

test_that("fit_initial() accepts model strings with or without whitespace", {
  compact <- list(
    z = list(
      model = "Surv(tstart,tstop,delta==2)~L0",
      at_risk = function(dt) dt[["z"]] == 0
    ),
    death = "Surv(tstart,tstop,delta==1)~L0+z",
    censoring = "Surv(tstart,tstop,delta==0)~L0"
  )

  expect_equal(
    fit_initial(illness_death, fit_types = compact, tau = 1)$tmp_long,
    fit_initial(
      illness_death,
      fit_types = illness.death.fit.types,
      tau = 1
    )$tmp_long
  )
})

test_that("fit_initial() accepts transformed covariates in models", {
  transformed <- illness.death.fit.types
  transformed$death <- "Surv(tstart, tstop, delta == 1) ~ log(exp(L0)) + z"

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
      fit_types = c(ft, list(other = "Surv(tstart, tstop, event) ~ L0")),
      tau = 1
    ),
    "delta == <code>"
  )
  expect_error(
    fit_initial(
      illness_death,
      fit_types = list(
        z = list(model = "Surv(tstart, tstop, delta == 2) ~ L0", fit = "glm")
      ),
      tau = 1
    ),
    "fit"
  )
  expect_error(
    fit_initial(illness_death, fit_types = ft, tau = 1, a = 1),
    "treatment_model"
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
      treatment_model = list(model = "A0 ~ L0", fit = "other")
    ),
    "glm"
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
    treatment_model = "A0 ~ L0",
    a = 0:1
  )

  expect_true(data.table::is.data.table(fit$tmp_long))
  expect_true(all(fit$tmp_long$pi.A0.1 >= 0 & fit$tmp_long$pi.A0.1 <= 1))
  expect_true(all(
    c("clever.weight.a0", "clever.weight.a1") %in%
      names(fit$tmp_long)
  ))
})

test_that("fit_initial() prints fit diagnostics when verbose = TRUE", {
  expect_output(suppressMessages(fixture.fit(20, verbose = TRUE)))
})

test_that("fit_initial()'s HAL setup identifies the process for its own delta", {
  # Regression test for a `name.hal` indexing bug that used to select every
  # process name instead of just the HAL one and error before reaching the
  # HAL fit itself.
  ft <- illness.death.fit.types
  ft$z$fit <- "hal"

  expect_error(
    fit_initial(fixture.data(30), fit_types = ft, tau = 1),
    "could not find function \"fit.hal\"",
    fixed = TRUE
  )
})
