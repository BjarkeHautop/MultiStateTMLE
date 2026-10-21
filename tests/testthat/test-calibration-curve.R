test_that("calibration_curve() returns one row per alpha", {
  fit <- fit_initial(
    illness_death,
    fit_types = illness.death.fit.types,
    tau = 1
  )

  curve <- calibration_curve(
    fit,
    target = "death",
    alpha_grid = c(0, 1),
    output_eic = TRUE
  )

  expect_equal(curve$estimate$alpha, c(0, 1))
  expect_named(
    curve$estimate,
    c("alpha", "exposure.est", "exposure.se", "target.est", "target.se")
  )
  # Without any z events, the exposure risk is zero.
  expect_equal(curve$estimate$exposure.est[1], 0)
  # alpha = 1 is no intervention.
  expect_equal(
    curve$estimate$target.est[2],
    tmle_alpha(fit, target = "death", min_iter = 1)$estimate[["tmle.est"]]
  )
  expect_named(curve$eics, c("alpha = 0", "alpha = 1"))
  expect_false(anyNA(curve$eics[["alpha = 0"]]$eic))
})

test_that("calibration_curve() validates its arguments", {
  fit <- fixture.fit(30)

  expect_error(
    calibration_curve(fit, target = "death", alpha_grid = numeric(0)),
    "alpha_grid"
  )
  expect_error(calibration_curve(fit, target = 1), "target")
  expect_error(
    calibration_curve(fit, target = "death", min_iter = -1),
    "min_iter"
  )
})

test_that("calibration_curve() reports progress only when verbose = TRUE", {
  fit <- fixture.fit(30)

  msgs <- capture_messages(
    calibration_curve(fit, target = "death", alpha_grid = 1, verbose = TRUE)
  )
  expect_true(any(grepl("TMLE iteration 1", msgs)))
  expect_silent(
    calibration_curve(fit, target = "death", alpha_grid = 1)
  )
})
