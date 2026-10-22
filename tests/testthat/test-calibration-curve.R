test_that("calibration_curve() returns one row per alpha", {
  fit <- fit_initial(
    illness_death,
    fit_types = illness.death.fit.types,
    tau = 1
  )

  curve <- calibration_curve(
    fit,
    z_name = "z",
    target = c("death", "z"),
    alpha_grid = c(0, 1)
  )

  expect_s3_class(curve, "calibration_curve")
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
    tmle_alpha(fit, target = "death")$estimate$estimate
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

test_that("calibration_curve() validates its arguments", {
  fit <- fixture.fit(30)

  expect_error(
    calibration_curve(
      fit,
      z_name = "z",
      target = "death",
      alpha_grid = numeric(0)
    ),
    "alpha_grid"
  )
  expect_error(calibration_curve(fit, z_name = "z", target = 1), "target")
  expect_error(
    calibration_curve(
      fit,
      z_name = "z",
      target = "death",
      control = list(min_iter = -1)
    ),
    "control\\$min_iter"
  )
  expect_error(
    calibration_curve(fit, z_name = "z", target = "death", min_iter = 1),
    "unused argument"
  )
})

test_that("calibration_curve() reports progress only when verbose = TRUE", {
  fit <- fixture.fit(30)

  msgs <- capture_messages(
    calibration_curve(
      fit,
      z_name = "z",
      target = "death",
      alpha_grid = 1,
      verbose = TRUE
    )
  )
  expect_true(any(grepl("TMLE iteration 1", msgs)))
  expect_silent(
    calibration_curve(fit, z_name = "z", target = "death", alpha_grid = 1)
  )
})

test_that("calibration_curve() has separate exposure and target estimands", {
  fit <- fixture.fit(30)

  curve <- calibration_curve(
    fit,
    z_name = "z",
    target = "death",
    alpha_grid = 1,
    estimand = "years_lost"
  )

  expect_equal(
    curve$estimate$exposure_estimate,
    tmle_alpha(fit, target = "z")$estimate$estimate
  )
  expect_equal(
    curve$estimate$estimate,
    tmle_alpha(
      fit,
      target = "death",
      estimand = "years_lost"
    )$estimate$estimate
  )
})
