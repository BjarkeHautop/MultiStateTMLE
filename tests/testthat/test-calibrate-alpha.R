test_that("calibrate_alpha() hits the calibration target", {
  fit <- fit_initial(
    illness_death,
    fit_types = illness.death.fit.types,
    tau = 1
  )
  n <- length(unique(illness_death$id))

  cal <- calibrate_alpha(fit, target = "outcome1", rho = 0.8, output_eic = TRUE)

  expect_named(
    cal$estimate,
    c(
      "alpha.est",
      "alpha.se",
      "target.est",
      "target.se",
      "target.se.crude",
      "theta",
      "theta.se"
    )
  )
  # Decreasing the risk of z requires decreasing its intensity.
  expect_lt(cal$estimate[["alpha.est"]], 1)

  # solve_alpha() stops once within c_n = n^(-1/2) / log(n) of theta.
  z_risk <- tmle_alpha(fit, target = "z", alpha = cal$estimate[["alpha.est"]])
  expect_lt(
    abs(z_risk$estimate[["tmle.est"]] - cal$estimate[["theta"]]),
    n^(-1 / 2) / log(n)
  )
  expect_length(cal$eic, n)
})

test_that("calibrate_alpha() validates its arguments", {
  fit <- fit_initial(
    illness_death,
    fit_types = illness.death.fit.types,
    tau = 1
  )

  expect_error(calibrate_alpha(fit, target = "outcome1"), "exactly one")
  expect_error(
    calibrate_alpha(fit, target = "outcome1", theta = 0.2, rho = 0.5),
    "exactly one"
  )
  expect_error(
    calibrate_alpha(fit, target = "outcome1", rho_1a = 0.5),
    "requires a baseline treatment"
  )
})
