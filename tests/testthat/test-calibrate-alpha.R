test_that("calibrate_alpha() hits the calibration target", {
  fit <- fit_initial(
    illness_death,
    fit_types = illness.death.fit.types,
    tau = 1
  )
  n <- length(unique(illness_death$id))

  cal <- calibrate_alpha(fit, z_name = "z", target = "death", rho = 0.8)

  expect_s3_class(cal, "calibrate_alpha")
  expect_named(
    cal$estimate,
    c("target", "estimate", "se", "lower", "upper", "se_crude")
  )
  expect_equal(cal$calibration$quantity, c("alpha", "theta"))
  alpha_hat <- cal$calibration$estimate[1]
  theta <- cal$calibration$estimate[2]
  # Decreasing the risk of z requires decreasing its intensity.
  expect_lt(alpha_hat, 1)

  # solve_alpha() stops once within c_n = n^(-1/2) / log(n) of theta.
  z_risk <- tmle_alpha(fit, target = "z", alpha = alpha_hat, z_name = "z")
  expect_lt(abs(z_risk$estimate$estimate - theta), n^(-1 / 2) / log(n))
  expect_named(cal$eic, "death")
  expect_length(cal$eic$death, n)
  expect_length(cal$alpha_eic, n)
})

test_that("calibrate_alpha() validates its arguments", {
  fit <- fit_initial(
    illness_death,
    fit_types = illness.death.fit.types,
    tau = 1
  )

  expect_error(
    calibrate_alpha(fit, z_name = "z", target = "death"),
    "exactly one"
  )
  expect_error(
    calibrate_alpha(
      fit,
      z_name = "z",
      target = "death",
      theta = 0.2,
      rho = 0.5
    ),
    "exactly one"
  )
  expect_error(
    calibrate_alpha(fit, z_name = "z", target = "death", rho_1a = 0.5),
    "requires a baseline treatment"
  )
})

test_that("calibrate_alpha() calibrates to theta directly", {
  fit <- fixture.fit(100)

  cal <- calibrate_alpha(fit, z_name = "z", target = "death", theta = 0.2)

  expect_equal(cal$calibration$estimate[2], 0.2)
  expect_true(cal$checks$converged)
  z_risk <- tmle_alpha(
    fit,
    target = "z",
    alpha = cal$calibration$estimate[1],
    z_name = "z"
  )
  expect_lt(abs(z_risk$estimate$estimate - 0.2), 0.02)
})

test_that("calibrate_alpha() handles several targets", {
  fit <- fixture.fit(100)

  single <- calibrate_alpha(fit, z_name = "z", target = "death", theta = 0.2)
  both <- calibrate_alpha(
    fit,
    z_name = "z",
    target = c("death", "z"),
    theta = 0.2
  )

  expect_equal(both$estimate$target, c("death", "z"))
  expect_equal(both$estimate[1], single$estimate)
  expect_equal(both$estimate$estimate[2], 0.2, tolerance = 0.1)
  expect_named(both$eic, c("death", "z"))
  expect_named(both$checks$derivative_target, c("death", "z"))
  expect_equal(coef(both), c(death = single$estimate$estimate, z = 0.2),
    tolerance = 0.1
  )
})

test_that("calibrate_alpha() rejects invalid inputs", {
  fit <- fixture.fit(30)

  expect_error(
    calibrate_alpha(fit, z_name = "z", target = "death", theta = "x"),
    "theta"
  )
  expect_error(
    calibrate_alpha(fit, z_name = "z", target = 5, theta = 0.2),
    "target"
  )
  expect_error(
    calibrate_alpha(
      fit,
      z_name = "z",
      target = "death",
      theta = 0.2,
      cores = 0
    ),
    "cores"
  )
})
