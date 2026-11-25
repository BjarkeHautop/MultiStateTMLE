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
      name = c("illness", "death", "censoring"),
      type = c("one_jump", "terminal", "terminal"),
      model = c(
        "illness ~ biomarker",
        "death ~ biomarker + illness",
        "censoring ~ biomarker"
      ),
      at_risk = c(TRUE, FALSE, FALSE)
    )
  )
  expect_named(fit$models, c("illness", "death", "censoring"))
  expect_s3_class(fit$models$death, "hazard_model")

  expect_equal(fit$process_names, c("illness", "death", "censoring"))
  expect_equal(
    fit$process_types,
    list(illness = "one_jump", death = "terminal", censoring = "terminal")
  )
  expect_equal(fit$cens_process_id, 3)
  expect_named(fit$at_risks, "illness")
  # illness is a state variable, with states illness = 0 and illness = 1.
  expect_equal(fit$depend_matrix$illness, c(0, 1))
  expect_true(all(fit$tmp_long$time <= 1))
  expect_true(all(
    c("P.illness.1", "P.illness.2", "P.death.1", "P.death.2") %in%
      names(fit$tmp_long)
  ))
})

test_that("fit_initial() takes simevent output as is", {
  skip_if_not_installed("simevent")
  sim <- simulate.illness.death(60, seed = 1)

  # simevent's illness column, the count after each event, is recomputed.
  expect_true("illness" %in% names(sim))
  fit <- fit_initial(sim, models = illness.death.models, tau = 1)
  expect_equal(fit$n, 60)
  expect_equal(fit$depend_matrix$illness, c(0, 1))
})

test_that("fit_initial() accepts plain formulas, named lists and any order", {
  fit <- fixture.fit(40)

  reordered <- fit_initial(
    fixture.data(40),
    models = list(
      censoring = censoring ~ biomarker,
      death = death ~ biomarker + illness,
      illness = hazard_model(illness ~ biomarker, at_risk = ~ illness == 0)
    ),
    tau = 1
  )
  expect_equal(
    tmle_estimate(
      reordered,
      risk(death),
      scale_intensity(illness = 0.5)
    )$estimate,
    tmle_estimate(fit, risk(death), scale_intensity(illness = 0.5))$estimate
  )
})

test_that("fit_initial() accepts transformed covariates in models", {
  transformed <- illness.death.models
  transformed[[2]] <- death ~ log(exp(biomarker)) + illness

  expect_equal(
    fit_initial(fixture.data(40), models = transformed, tau = 1)$tmp_long,
    fixture.fit(40)$tmp_long
  )
})

test_that("fit_initial() derives covariates from processes with I()", {
  # With illness happening at most once, I(illness > 0) equals illness.
  derived <- illness.death.models
  derived[[2]] <- death ~ biomarker + I(illness > 0)
  fit <- fit_initial(fixture.data(40), models = derived, tau = 1)

  expect_equal(fit$depend_matrix$illness, c(0, 1))
  expect_equal(
    tmle_estimate(fit, risk(death), scale_intensity(illness = 0.5))$estimate,
    tmle_estimate(
      fixture.fit(40),
      risk(death),
      scale_intensity(illness = 0.5)
    )$estimate
  )
})

test_that("fit_initial() uses a censoring process of another name", {
  data <- fixture.data(40)
  levels(data$event)[levels(data$event) == "censoring"] <- "dropout"

  expect_error(
    fit_initial(
      data,
      models = list(
        hazard_model(illness ~ biomarker, limit = 1),
        death ~ biomarker + illness
      ),
      tau = 1
    ),
    "censoring"
  )
  fit <- fit_initial(
    data,
    models = list(
      hazard_model(illness ~ biomarker, limit = 1),
      death ~ biomarker + illness,
      dropout ~ biomarker
    ),
    tau = 1,
    censoring = "dropout"
  )
  expect_equal(fit$process_names[fit$cens_process_id], "dropout")
})

test_that("hazard_model() validates its arguments", {
  expect_error(hazard_model(~biomarker), "two-sided")
  expect_error(hazard_model(log(illness) ~ biomarker), "two-sided")
  expect_error(hazard_model(illness ~ biomarker, limit = 0), "limit")
  expect_error(
    hazard_model(illness ~ biomarker, at_risk = illness ~ illness == 0),
    "one-sided"
  )
  expect_error(hazard_model(illness ~ biomarker, fit = "glm"), "cox")
  expect_output(
    print(hazard_model(illness ~ biomarker, limit = 1, at_risk = ~ death == 0)),
    "limit: 1"
  )
})

test_that("hazard_model() conditions become at-risk functions", {
  k <- 5
  at_risk <- .hazard_model_to_fit_type(
    hazard_model(
      illness ~ biomarker,
      limit = 1,
      at_risk = ~ death == 0 & illness < k
    ),
    delta = 1,
    state_vars = c("illness", "death")
  )$at_risk

  # State variables are read from dt, other variables from the formula's
  # environment.
  expect_equal(
    all.names(body(at_risk)),
    c("&", "<", "[[", "dt", "&", "==", "[[", "dt", "<", "[[", "dt", "k")
  )
  dt <- data.table::data.table(illness = c(0, 1, 0), death = c(0, 0, 1))
  expect_equal(at_risk(dt), c(TRUE, FALSE, FALSE))
  expect_null(
    .hazard_model_to_fit_type(
      hazard_model(illness ~ biomarker),
      1,
      "illness"
    )$at_risk
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
    fit_initial(illness_death, models = c(m, list(death ~ biomarker)), tau = 1),
    "one model"
  )
  expect_error(
    fit_initial(illness_death, models = c(m, list("x")), tau = 1),
    "formula or a hazard_model"
  )
  expect_error(
    fit_initial(
      illness_death,
      models = list(other = death ~ biomarker),
      tau = 1
    ),
    "named \"other\""
  )
  expect_error(
    fit_initial(
      illness_death,
      models = list(
        hazard_model(illness ~ biomarker, at_risk = ~ biomarker > 0),
        death ~ biomarker + illness,
        censoring ~ biomarker
      ),
      tau = 1
    ),
    "biomarker"
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
      fixture.data(),
      models = m,
      tau = 1,
      treatment_model = list(model = treatment ~ biomarker, fit = "glm")
    ),
    "treatment_model"
  )
  expect_error(
    fit_initial(
      fixture.data(),
      models = m,
      tau = 1,
      treatment_model = ~biomarker,
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
    fixture.data(),
    models = illness.death.models,
    tau = 1,
    treatment_model = treatment ~ biomarker
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
  m[[1]] <- hazard_model(illness ~ biomarker, limit = 1, fit = "hal")
  assign("fit.hal", function(...) stop("reached fit.hal"), envir = globalenv())
  withr::defer(rm("fit.hal", envir = globalenv()))

  expect_error(
    fit_initial(fixture.data(30), models = m, tau = 1),
    "reached fit.hal"
  )
})

test_that("fit_initial() errors early when HAL is not available", {
  m <- illness.death.models
  m[[1]] <- hazard_model(illness ~ biomarker, limit = 1, fit = "hal")

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

test_that("fit_initial() takes event_time() terms in models", {
  data <- fixture.data(40)
  timed <- fit_initial(
    data,
    models = list(
      hazard_model(illness ~ biomarker, limit = 1),
      death ~ biomarker + illness + I(event_time(illness, 1) > 0.5),
      censoring ~ biomarker + I(event_time(illness) > 0.5)
    ),
    tau = 1
  )
  expect_equal(
    timed$processes$model[2],
    "death ~ biomarker + illness + I(event_time(illness, 1) > 0.5)"
  )

  # The term is 1 from a first illness after time 0.5 on, and 0 otherwise.
  data.table::setDT(data)
  first <- data[event == "illness", list(first = min(time)), by = "id"]
  rows <- merge(timed$tmp_long, first, by = "id", all.x = TRUE)
  expected <- !is.na(rows$first) & rows$first > 0.5 & rows$tstart >= rows$first
  expect_equal(rows$T.illness, as.numeric(expected))

  # Categories of the time, with one effect each.
  categories <- fit_initial(
    data,
    models = list(
      hazard_model(illness ~ biomarker, limit = 1),
      death ~ biomarker +
        illness +
        I((event_time(illness) > 0) + (event_time(illness) > 0.5)),
      censoring ~ biomarker
    ),
    tau = 1
  )
  expect_setequal(unique(categories$tmp_long$T.illness), 0:2)
  expect_length(coef(categories$cox_fits$death), 4)

  # The process of the event time is a state variable even when its count
  # is not in any model.
  no_count <- fit_initial(
    data,
    models = list(
      hazard_model(illness ~ biomarker, limit = 1),
      death ~ biomarker + I(event_time(illness) > 0.5),
      censoring ~ biomarker
    ),
    tau = 1
  )
  expect_true("illness" %in% names(no_count$depend_matrix))
})

test_that("fit_initial() takes event_time() terms that are 1 before the event", {
  fit <- function(term) {
    fit_initial(
      fixture.data(40),
      models = list(
        hazard_model(illness ~ biomarker, limit = 1),
        stats::as.formula(paste("death ~ biomarker + illness +", term)),
        censoring ~ biomarker
      ),
      tau = 1
    )
  }
  # Complementary terms: the same model, with the sign of the effect flipped.
  early <- fit("I(event_time(illness) < 0.5)")
  late <- fit("I(event_time(illness) > 0 & event_time(illness) >= 0.5)")

  expect_equal(
    coef(early$cox_fits$death)[["as.factor(T.illness)1"]],
    -coef(late$cox_fits$death)[["as.factor(T.illness)1"]]
  )
  expect_equal(early$tmp_long$P.death, late$tmp_long$P.death)
  # Before the illness, the term is 1.
  expect_true(all(early$depend_matrix[illness == 0]$T.illness == 1))
  expect_equal(
    tmle_estimate(early, risk(death), scale_intensity(illness = 0.5))$estimate,
    tmle_estimate(late, risk(death), scale_intensity(illness = 0.5))$estimate
  )
})

test_that("fit_initial() requires event_time() terms to be categories", {
  models <- function(term) {
    list(
      hazard_model(illness ~ biomarker, limit = 1),
      stats::as.formula(paste("death ~ biomarker +", term)),
      censoring ~ biomarker
    )
  }
  fit <- function(term) {
    fit_initial(fixture.data(40), models = models(term), tau = 1)
  }
  expect_error(fit("event_time(illness)"), "is the event time itself")
  expect_error(
    fit("I(exp(-event_time(illness)))"),
    "values other than whole numbers"
  )
  expect_error(fit("I(event_time(illness, 2) > 0.5)"), "max_count")
  expect_error(fit("I(event_time(foo) > 0.5)"), "not a process")
  expect_error(fit("I(event_time(c(illness, foo)) > 0.5)"), "foo, which is not")
  expect_error(fit("I(event_time('illness') > 0.5)"), "give a process")
  expect_error(fit("I(event_time(illness) > biomarker)"), "only use constants")
  expect_error(
    fit("I(event_time(illness) > 0.5) + I(event_time(illness) > 0.7)"),
    "several functions"
  )
})

# Two processes, a and b, that each happen at most once, b only after a,
# and death, whose hazard is lower when a happened late.
simulate.two.jumps <- function(n, seed) {
  model <- simevent::sim_model(
    biomarker = simevent::sim_covariate(function(N) rnorm(N)),
    a = simevent::sim_process("transient", limit = 1, eta = 0.6, nu = 1),
    b = simevent::sim_process("transient", limit = 1, eta = 0.8, nu = 1),
    death = simevent::sim_process("terminal", eta = 0.3, nu = 1),
    censoring = simevent::sim_process("censoring", eta = 0.2, nu = 1),
    effects = list(
      simevent::sim_effect(b ~ (a == 0), coef = -50),
      simevent::sim_effect(death ~ biomarker + a, coef = c(0.5, 0.7))
    )
  )
  simevent::simulate_events(model, n = n, max_cens = 2, seed = seed)
}

two.jumps.models <- function(term) {
  list(
    hazard_model(a ~ biomarker, limit = 1),
    hazard_model(b ~ biomarker, limit = 1, at_risk = ~ a == 1),
    stats::as.formula(paste("death ~ biomarker + a +", term)),
    censoring ~ biomarker
  )
}

test_that("event_time() takes several processes", {
  skip_if_not_installed("simevent")
  data <- data.table::as.data.table(simulate.two.jumps(150, seed = 1))
  data[event == "max_cens", event := "censoring"]
  fit <- function(term) {
    # max_count = 2 to track the second event of a or b.
    fit_initial(data, models = two.jumps.models(term), tau = 1.5, max_count = 2)
  }

  # The term is 1 from a first event of a or b after time 0.5 on.
  either <- fit("I(event_time(c(b, a)) > 0.5)")
  first <- data[event %in% c("a", "b"), list(first = min(time)), by = "id"]
  rows <- merge(either$tmp_long, first, by = "id", all.x = TRUE)
  expected <- !is.na(rows$first) & rows$first > 0.5 & rows$tstart >= rows$first
  expect_equal(rows$T.a_b, as.numeric(expected))
  expect_true(all(either$tmp_long$state > 0))

  # As b only happens after a, the second event of a or b is the first b,
  # and the first event of a or b is the first a: the same fits.
  second <- fit("I(event_time(c(a, b), 2) > 0.5)")
  only_b <- fit("I(event_time(b) > 0.5)")
  only_a <- fit("I(event_time(a) > 0.5)")
  expect_equal(
    unname(coef(second$cox_fits$death)),
    unname(coef(only_b$cox_fits$death))
  )
  expect_equal(
    unname(coef(either$cox_fits$death)),
    unname(coef(only_a$cox_fits$death))
  )
  estimate <- function(fit) {
    tmle_estimate(fit, risk(death), scale_intensity(a = 0.5))$estimate
  }
  expect_equal(estimate(second), estimate(only_b))
  expect_equal(estimate(either), estimate(only_a))
})

test_that("event_time() of several processes interacts with counts", {
  skip_if_not_installed("simevent")
  data <- data.table::as.data.table(simulate.two.jumps(150, seed = 2))
  data[event == "max_cens", event := "censoring"]
  fit <- fit_initial(
    data,
    models = two.jumps.models("I(event_time(c(a, b), 1) > 0.5) * I(b == 0)"),
    tau = 1.5
  )
  expect_named(
    coef(fit$cox_fits$death),
    c(
      "biomarker",
      "a",
      "as.factor(T.a_b)1",
      "I(b == 0)TRUE",
      "as.factor(T.a_b)1:I(b == 0)TRUE"
    )
  )
  estimate <- tmle_estimate(fit, risk(death))$estimate$estimate
  expect_true(estimate > 0 && estimate < 1)
})

test_that("fit_initial() takes factor and character covariates", {
  data <- fixture.data(40)
  numeric <- fit_initial(
    data,
    models = list(
      hazard_model(illness ~ biomarker + treatment, limit = 1),
      death ~ biomarker + treatment + illness,
      censoring ~ biomarker
    ),
    tau = 1
  )
  for (as_type in list(as.factor, as.character)) {
    typed <- data
    typed$treatment <- as_type(ifelse(data$treatment == 1, "yes", "no"))
    fit <- fit_initial(
      typed,
      models = list(
        hazard_model(illness ~ biomarker + treatment, limit = 1),
        death ~ biomarker + treatment + illness,
        censoring ~ biomarker
      ),
      tau = 1
    )
    expect_equal(
      unname(coef(fit$cox_fits$death)),
      unname(coef(numeric$cox_fits$death))
    )
    expect_equal(fit$tmp_long$P.death, numeric$tmp_long$P.death)
  }
})
