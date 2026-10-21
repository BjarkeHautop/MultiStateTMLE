test_that("tmle_alpha() returns the estimate and requested extras", {
  fit <- fit_initial(
    illness_death,
    fit_types = illness.death.fit.types,
    tau = 1
  )

  est <- tmle_alpha(
    fit,
    target = "death",
    alpha = 0.5,
    output_eic = TRUE,
    output_convergence = TRUE,
    output_weights = c(0.5, 1)
  )

  expect_named(est$estimate, c("tmle_est", "se", "g_est"))
  expect_true(est$estimate[["tmle_est"]] > 0 && est$estimate[["tmle_est"]] < 1)
  expect_length(est$eic, length(unique(illness_death$id)))
  expect_named(est$convergence, c("iter", "eic_solved_at", "converged"))
  expect_named(est$weights, c("q50", "q100"))
  expect_named(est$cens_weights, c("q50", "q100"))

  one_step <- tmle_alpha(fit, target = "death", alpha = 0.5, one_step = TRUE)
  expect_named(one_step$estimate, c("one_step_est", "se", "g_est"))
  expect_equal(one_step$estimate[["g_est"]], est$estimate[["g_est"]])
})

test_that("tmle_alpha() matches the target case-insensitively", {
  fit <- fit_initial(
    illness_death,
    fit_types = illness.death.fit.types,
    tau = 1
  )

  expect_equal(
    tmle_alpha(fit, target = "DEATH", alpha = 0.5)$estimate,
    tmle_alpha(fit, target = "death", alpha = 0.5)$estimate
  )
})

test_that("tmle_alpha() validates its arguments", {
  fit <- fit_initial(
    illness_death,
    fit_types = illness.death.fit.types,
    tau = 1
  )

  expect_error(tmle_alpha(list(), target = "death"), "initial_fit")
  expect_error(tmle_alpha(fit, target = "stroke"), "target must be one of")
  expect_error(tmle_alpha(fit, target = "death", alpha = -1), "alpha")
  expect_error(tmle_alpha(fit, target = "death", a = 1), "clever.weight.a1")
  expect_warning(
    tmle_alpha(fit, target = "death", alpha = 0.5, z_name = "w"),
    "alpha is ignored"
  )
})

test_that("tmle_alpha() supports only_first, target_only_in_state, target_by_state", {
  fit <- fixture.fit(30)

  expect_true(is.numeric(
    tmle_alpha(fit, target = "z", only_first = TRUE)$estimate[["tmle_est"]]
  ))
  expect_true(is.numeric(
    tmle_alpha(
      fit,
      target = "death",
      target_only_in_state = function(states) states[["z"]] == 0
    )$estimate[["tmle_est"]]
  ))
  expect_true(is.numeric(
    tmle_alpha(fit, target = "death", target_by_state = TRUE)$estimate[[
      "tmle_est"
    ]]
  ))
})

test_that("tmle_alpha() truncates clever weights", {
  est <- tmle_alpha(
    fixture.fit(30),
    target = "death",
    alpha = 3,
    truncate_weights = 1.5,
    output_weights = 0.5
  )

  expect_lte(est$weights[["q50"]], 1.5 + 1e-8)
  expect_true("n_truncated" %in% names(est$weights))
})

test_that("tmle_alpha() estimates under a baseline intervention 'a'", {
  fit <- fit_initial(
    fixture.data.treated(),
    fit_types = illness.death.fit.types,
    tau = 1,
    treatment_model = "A0 ~ L0",
    a = 0:1
  )

  est <- tmle_alpha(fit, target = "death", a = 1, output_a_weights = 0.5)
  expect_named(est$a_weights, "q50")
  expect_true(
    est$estimate[["tmle_est"]] >= 0 && est$estimate[["tmle_est"]] <= 1
  )

  expect_error(tmle_alpha(fit, target = "death", a = 5), "a")
})

test_that("tmle_alpha()'s alpha_list reproduces a constant scalar alpha", {
  fit <- fixture.fit(30)
  constant <- list(z = function(time) rep(2, length(time)))

  scalar <- tmle_alpha(fit, target = "death", alpha = 2)
  expect_equal(
    tmle_alpha(fit, target = "death", alpha_list = constant)$estimate,
    scalar$estimate
  )
  expect_equal(
    tmle_alpha(
      fit,
      target = "death",
      alpha_list = constant,
      target_by_state = TRUE
    )$estimate[["tmle_est"]],
    scalar$estimate[["tmle_est"]]
  )
})
