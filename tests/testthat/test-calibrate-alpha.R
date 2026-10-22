test_that("calibrate_alpha() hits the calibration target", {
  fit <- fit_initial(
    illness_death,
    fit_types = illness.death.fit.types,
    tau = 1
  )
  n <- length(unique(illness_death$id))

  cal <- calibrate_alpha(
    fit,
    z_name = "z",
    target = "death",
    value = 0.8,
    scale = "ratio"
  )

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
    "value"
  )
  expect_error(
    calibrate_alpha(
      fit,
      z_name = "z",
      target = "death",
      value = 0.5,
      scale = "fraction"
    ),
    "should be one of"
  )
  expect_error(
    calibrate_alpha(
      fit,
      z_name = "z",
      target = "death",
      value = 0.5,
      scale = "ratio",
      reference_a = 1
    ),
    "reference_a"
  )
  expect_error(
    calibrate_alpha(
      fit,
      z_name = "z",
      target = "death",
      value = 0.5,
      estimand = "rate"
    ),
    "should be one of"
  )
})

test_that("calibrate_alpha() calibrates to theta directly", {
  fit <- fixture.fit(100)

  cal <- calibrate_alpha(fit, z_name = "z", target = "death", value = 0.2)

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

  single <- calibrate_alpha(fit, z_name = "z", target = "death", value = 0.2)
  both <- calibrate_alpha(
    fit,
    z_name = "z",
    target = c("death", "z"),
    value = 0.2
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
    calibrate_alpha(fit, z_name = "z", target = "death", value = "x"),
    "value"
  )
  expect_error(
    calibrate_alpha(fit, z_name = "z", target = 5, value = 0.2),
    "target"
  )
  expect_error(
    calibrate_alpha(
      fit,
      z_name = "z",
      target = "death",
      value = 0.2,
      cores = 0
    ),
    "cores"
  )
})

test_that("calibrate_alpha() calibrates relative to a reference", {
  fit <- fixture.fit(100)
  reference <- tmle_alpha(fit, target = "z")$estimate$estimate

  cal <- calibrate_alpha(
    fit,
    z_name = "z",
    target = "death",
    value = -0.05,
    scale = "difference"
  )
  expect_equal(cal$calibration$estimate[2], reference - 0.05)
  expect_equal(cal$scale, "difference")
})

test_that("calibrate_alpha() calibrates to the risk under another arm", {
  fit <- fit_initial(
    fixture.data.treated(100),
    fit_types = illness.death.fit.types,
    tau = 1,
    treatment_model = A0 ~ L0,
    a = 0:1
  )

  cal <- calibrate_alpha(
    fit,
    z_name = "z",
    target = "death",
    value = 1,
    scale = "ratio",
    a = 1,
    reference_a = 0
  )
  expect_equal(
    cal$calibration$estimate[2],
    tmle_alpha(fit, target = "z", a = 0)$estimate$estimate
  )
})

test_that("calibrate_alpha() can calibrate on risk and target years lost", {
  fit <- fixture.fit(100)

  cal <- calibrate_alpha(
    fit,
    z_name = "z",
    target = "death",
    value = 0.8,
    scale = "ratio",
    estimand = "years_lost"
  )

  # The exposure is calibrated on the risk scale...
  expect_equal(
    cal$calibration$estimate[2],
    0.8 * tmle_alpha(fit, target = "z")$estimate$estimate
  )
  # ...while the target is on the years-lost scale.
  expect_equal(
    cal$estimate$estimate,
    tmle_alpha(
      fit,
      target = "death",
      alpha = cal$calibration$estimate[1],
      z_name = "z",
      estimand = "years_lost"
    )$estimate$estimate
  )
  expect_equal(cal$estimand, "years_lost")
  expect_equal(cal$exposure_estimand, "risk")
})
