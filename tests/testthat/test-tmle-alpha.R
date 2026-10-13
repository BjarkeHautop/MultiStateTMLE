test_that("tmle_alpha() returns the estimate and requested extras", {
  fit <- fit_initial(
    illness_death,
    fit_types = illness.death.fit.types,
    tau = 1
  )

  est <- tmle_alpha(
    fit,
    target = "outcome1",
    alpha = 0.5,
    output_eic = TRUE,
    output_convergence = TRUE,
    output_weights = c(0.5, 1)
  )

  expect_named(est$estimate, c("tmle.est", "se", "g.est"))
  expect_true(est$estimate[["tmle.est"]] > 0 && est$estimate[["tmle.est"]] < 1)
  expect_length(est$eic, length(unique(illness_death$id)))
  expect_named(est$convergence, c("iter", "eic.solved.at", "converged"))
  expect_named(est$weights, c("q50", "q100"))
  expect_named(est$cens.weights, c("q50", "q100"))

  one_step <- tmle_alpha(fit, target = "outcome1", alpha = 0.5, one_step = TRUE)
  expect_named(one_step$estimate, c("one.step.est", "se", "g.est"))
  expect_equal(one_step$estimate[["g.est"]], est$estimate[["g.est"]])
})

test_that("tmle_alpha() matches the target case-insensitively", {
  fit <- fit_initial(
    illness_death,
    fit_types = illness.death.fit.types,
    tau = 1
  )

  expect_equal(
    tmle_alpha(fit, target = "OUTCOME1", alpha = 0.5)$estimate,
    tmle_alpha(fit, target = "outcome1", alpha = 0.5)$estimate
  )
})

test_that("tmle_alpha() validates its arguments", {
  fit <- fit_initial(
    illness_death,
    fit_types = illness.death.fit.types,
    tau = 1
  )

  expect_error(tmle_alpha(list(), target = "outcome1"), "initial_fit")
  expect_error(tmle_alpha(fit, target = "death"), "target must be one of")
  expect_error(tmle_alpha(fit, target = "outcome1", alpha = -1), "alpha")
  expect_error(tmle_alpha(fit, target = "outcome1", a = 1), "clever.weight.a1")
  expect_warning(
    tmle_alpha(fit, target = "outcome1", alpha = 0.5, z_name = "w"),
    "alpha is ignored"
  )
})
