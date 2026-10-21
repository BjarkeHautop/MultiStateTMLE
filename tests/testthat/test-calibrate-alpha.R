test_that("calibrate_alpha() hits the calibration target", {
  fit <- fit_initial(
    illness_death,
    fit_types = illness.death.fit.types,
    tau = 1
  )
  n <- length(unique(illness_death$id))

  cal <- calibrate_alpha(fit, target = "death", rho = 0.8, output_eic = TRUE)

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
  expect_equal(cal$checks[["est.alpha.converged"]], 1)
  z_risk <- tmle_alpha(fit, target = "z", alpha = cal$estimate[["alpha.est"]])
  expect_lt(abs(z_risk$estimate[["tmle.est"]] - 0.2), 0.02)
})

test_that("calibrate_alpha() handles several targets", {
  fit <- fixture.fit(100)

  single <- calibrate_alpha(fit, target = "death", theta = 0.2)
  both <- calibrate_alpha(
    fit,
    target = c("death", "z"),
    theta = 0.2,
    output_eic = TRUE
  )

  expect_equal(
    both$estimate[["target.est.death"]],
    single$estimate[["target.est"]]
  )
  expect_equal(both$estimate[["target.est.z"]], 0.2, tolerance = 0.1)
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
