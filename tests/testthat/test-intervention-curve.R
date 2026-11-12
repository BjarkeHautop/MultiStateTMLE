test_that("intervention_curve() returns one row per alpha", {
  fit <- fit_initial(
    illness_death,
    models = illness.death.models,
    tau = 1
  )

  curve <- intervention_curve(fit, risk(death, z), risk(z), alpha = c(0, 1))

  expect_s3_class(curve, "intervention_curve")
  expect_equal(curve$estimate$alpha, c(0, 0, 1, 1))
  expect_equal(curve$estimate$target, c("death", "z", "death", "z"))
  expect_named(
    curve$estimate,
    c(
      "alpha",
      "target",
      "exposure_estimate",
      "exposure_se",
      "estimate",
      "se",
      "lower",
      "upper"
    )
  )
  # Without any z events, the exposure risk is zero.
  expect_equal(curve$estimate$exposure_estimate[1:2], c(0, 0))
  # alpha = 1 is no intervention.
  expect_equal(
    curve$estimate$estimate[3],
    tmle_estimate(fit, risk(death))$estimate$estimate
  )
  # The exposure is itself a target here.
  expect_equal(curve$estimate$estimate[4], curve$estimate$exposure_estimate[4])
  expect_equal(nrow(curve$checks), 4)
  expect_named(curve$eic, c("alpha = 0", "alpha = 1"))
  expect_named(
    curve$eic[["alpha = 1"]],
    c("exposure_eic", "alpha_eic", "eic", "eic_alpha_known")
  )
  expect_named(curve$eic[["alpha = 1"]]$eic, c("death", "z"))
  expect_false(anyNA(curve$eic[["alpha = 0"]]$eic$death))
})

test_that("intervention_curve() validates its arguments", {
  fit <- fixture.fit(30)

  expect_error(
    intervention_curve(fit, risk(death), risk(z), alpha = numeric(0)),
    "alpha"
  )
  expect_error(intervention_curve(fit, "death", risk(z)), "target")
  expect_error(intervention_curve(fit, risk(death), "z"), "exposure")
  expect_error(
    intervention_curve(fit, risk(death), risk(z, death)),
    "single process"
  )
  expect_error(
    intervention_curve(fit, risk(death), risk(w)),
    "w.*not in initial_fit"
  )
  expect_error(
    intervention_curve(
      fit,
      risk(death),
      risk(z),
      control = list(min_iter = -1)
    ),
    "control\\$min_iter"
  )
})

test_that("intervention_curve() reports progress only when verbose = TRUE", {
  fit <- fixture.fit(30)

  msgs <- capture_messages(
    intervention_curve(fit, risk(death), risk(z), alpha = 1, verbose = TRUE)
  )
  expect_true(any(grepl("TMLE iteration 1", msgs)))
  expect_silent(intervention_curve(fit, risk(death), risk(z), alpha = 1))
})

test_that("intervention_curve() has separate exposure and target estimands", {
  fit <- fixture.fit(30)

  curve <- intervention_curve(fit, years_lost(death), risk(z), alpha = 1)

  expect_equal(
    curve$estimate$exposure_estimate,
    tmle_estimate(fit, risk(z))$estimate$estimate
  )
  expect_equal(
    curve$estimate$estimate,
    tmle_estimate(fit, years_lost(death))$estimate$estimate
  )
})
