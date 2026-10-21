test_that("calibrate_alpha() hits the calibration target", {
  fit <- fit_initial(
    illness_death,
    fit_types = illness.death.fit.types,
    tau = 1
  )
  n <- length(unique(illness_death$id))

  cal <- calibrate_alpha(fit, target = "death", rho = 0.8)

  expect_named(
    cal$estimate,
    c(
      "alpha_est",
      "alpha_se",
      "target_est",
      "target_se",
      "target_se_crude",
      "theta",
      "theta_se"
    )
  )
  # Decreasing the risk of z requires decreasing its intensity.
  expect_lt(cal$estimate[["alpha_est"]], 1)

  # solve_alpha() stops once within c_n = n^(-1/2) / log(n) of theta.
  z_risk <- tmle_alpha(fit, target = "z", alpha = cal$estimate[["alpha_est"]])
  expect_lt(
    abs(z_risk$estimate[["tmle_est"]] - cal$estimate[["theta"]]),
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

  expect_error(calibrate_alpha(fit, target = "death"), "exactly one")
  expect_error(
    calibrate_alpha(fit, target = "death", theta = 0.2, rho = 0.5),
    "exactly one"
  )
  expect_error(
    calibrate_alpha(fit, target = "death", rho_1a = 0.5),
    "requires a baseline treatment"
  )
})

test_that("calibrate_alpha() calibrates to theta directly", {
  fit <- fixture.fit(100)

  cal <- calibrate_alpha(fit, target = "death", theta = 0.2)

  expect_equal(cal$estimate[["theta"]], 0.2)
  expect_equal(cal$checks[["est_alpha_converged"]], 1)
  z_risk <- tmle_alpha(fit, target = "z", alpha = cal$estimate[["alpha_est"]])
  expect_lt(abs(z_risk$estimate[["tmle_est"]] - 0.2), 0.02)
})

test_that("calibrate_alpha() handles several targets", {
  fit <- fixture.fit(100)

  single <- calibrate_alpha(fit, target = "death", theta = 0.2)
  both <- calibrate_alpha(fit, target = c("death", "z"), theta = 0.2)

  expect_equal(
    both$estimate[["target_est_death"]],
    single$estimate[["target_est"]]
  )
  expect_equal(both$estimate[["target_est_z"]], 0.2, tolerance = 0.1)
  expect_named(both$eic, c("death", "z"))
})

test_that("calibrate_alpha() rejects invalid inputs", {
  fit <- fixture.fit(30)

  expect_error(calibrate_alpha(fit, target = "death", theta = "x"), "theta")
  expect_error(calibrate_alpha(fit, target = 5, theta = 0.2), "target")
  expect_error(
    calibrate_alpha(fit, target = "death", theta = 0.2, cores = 0),
    "cores"
  )
})
