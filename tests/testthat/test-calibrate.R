test_that("calibrate() hits the calibration target", {
  fit <- fit_initial(
    illness_death,
    models = illness.death.models,
    tau = 1
  )
  n <- length(unique(illness_death$id))

  cal <- tmle_estimate(
    fit,
    risk(death),
    calibrate(risk(illness), relative_change(0.8))
  )

  expect_s3_class(cal, "tmle_estimate")
  expect_named(
    cal$estimate,
    c("target", "estimate", "se", "lower", "upper", "se_crude")
  )
  expect_equal(cal$calibration$quantity, c("alpha", "theta"))
  alpha_hat <- cal$calibration$estimate[1]
  theta <- cal$calibration$estimate[2]
  # Decreasing the risk of illness requires decreasing its intensity.
  expect_lt(alpha_hat, 1)

  # solve_alpha() stops once within c_n = n^(-1/2) / log(n) of theta.
  z_risk <- tmle_estimate(
    fit,
    risk(illness),
    scale_intensity(illness = alpha_hat)
  )
  expect_lt(abs(z_risk$estimate$estimate - theta), n^(-1 / 2) / log(n))
  expect_named(cal$eic, "death")
  expect_length(cal$eic$death, n)
  expect_length(cal$alpha_eic, n)
  expect_true(cal$checks$calibration$converged)
  expect_output(print(cal), "Estimates under the calibrated intervention")
})

test_that("calibrate() and the calibration rules validate their arguments", {
  expect_error(
    calibrate("illness", relative_change(0.8)),
    "risk\\(\\) or years_lost"
  )
  expect_error(calibrate(risk(illness, death), level(0.2)), "single process")
  expect_error(calibrate(risk(illness), 0.8), "level\\(\\), absolute_change")
  expect_error(
    calibrate(risk(illness), level(0.2), reference_a = 0),
    "reference_a"
  )
  expect_error(
    calibrate(risk(illness), level(0.2), reference_a = 2),
    "reference_a"
  )
  expect_error(level(1.5), "theta")
  expect_error(level("x"), "theta")
  expect_error(relative_change(-1), "rho")
  expect_error(absolute_change(Inf), "delta")
  expect_output(
    print(calibrate(risk(illness), absolute_change(-0.05), reference_a = 0)),
    "calibrate\\(risk\\(illness\\), absolute_change\\(-0.05\\), reference_a = 0\\)"
  )
})

test_that("calibrate() validates against the initial fit", {
  fit <- fixture.fit(30)

  expect_error(
    tmle_estimate(fit, risk(death), calibrate(risk(w), level(0.2))),
    "intervention is on w"
  )
  expect_error(
    tmle_estimate(
      fit,
      risk(death),
      calibrate(risk(illness), relative_change(0.8), reference_a = 1)
    ),
    "reference_a"
  )
  expect_error(
    tmle_estimate(
      fit,
      risk(death),
      calibrate(risk(illness), level(0.2)),
      control = list(max_iter = 0)
    ),
    "control\\$max_iter"
  )
})

test_that("calibrate() calibrates to a level", {
  fit <- fixture.fit(100)

  cal <- tmle_estimate(fit, risk(death), calibrate(risk(illness), level(0.2)))

  expect_equal(cal$calibration$estimate[2], 0.2)
  # A fixed level is known, so theta has no uncertainty.
  expect_equal(cal$calibration$se[2], 0)
  expect_true(cal$checks$calibration$converged)
  z_risk <- tmle_estimate(
    fit,
    risk(illness),
    scale_intensity(illness = cal$calibration$estimate[1])
  )
  expect_lt(abs(z_risk$estimate$estimate - 0.2), 0.02)
})

test_that("calibrate() handles several targets", {
  fit <- fixture.fit(100)
  rule <- calibrate(risk(illness), level(0.2))

  single <- tmle_estimate(fit, risk(death), rule)
  both <- tmle_estimate(fit, risk(death, illness), rule)

  expect_equal(both$estimate$target, c("death", "illness"))
  expect_equal(both$estimate[1], single$estimate)
  expect_equal(both$estimate$estimate[2], 0.2, tolerance = 0.1)
  expect_named(both$eic, c("death", "illness"))
  expect_named(both$checks$calibration$derivative_target, c("death", "illness"))
  expect_equal(
    coef(both),
    c(death = single$estimate$estimate, illness = 0.2),
    tolerance = 0.1
  )
})

test_that("calibrate() calibrates an absolute change", {
  fit <- fixture.fit(100)
  reference <- tmle_estimate(fit, risk(illness))$estimate$estimate

  cal <- tmle_estimate(
    fit,
    risk(death),
    calibrate(risk(illness), absolute_change(-0.05))
  )
  expect_equal(cal$calibration$estimate[2], reference - 0.05)
  expect_equal(cal$intervention$toward, absolute_change(-0.05))
})

test_that("calibrate() calibrates to the risk under another arm", {
  fit <- fit_initial(
    fixture.data(100),
    models = illness.death.models,
    tau = 1,
    treatment_model = treatment ~ biomarker,
    arms = 0:1
  )

  cal <- tmle_estimate(
    fit,
    risk(death),
    calibrate(risk(illness), relative_change(1), reference_a = 0),
    a = 1
  )
  expect_equal(
    cal$calibration$estimate[2],
    tmle_estimate(fit, risk(illness), a = 0)$estimate$estimate
  )
})

test_that("calibrate() can calibrate on risk and target years lost", {
  fit <- fixture.fit(100)

  cal <- tmle_estimate(
    fit,
    years_lost(death),
    calibrate(risk(illness), relative_change(0.8))
  )

  # The exposure is calibrated on the risk scale...
  expect_equal(
    cal$calibration$estimate[2],
    0.8 * tmle_estimate(fit, risk(illness))$estimate$estimate
  )
  # ...while the target is on the years-lost scale.
  expect_equal(
    cal$estimate$estimate,
    tmle_estimate(
      fit,
      years_lost(death),
      scale_intensity(illness = cal$calibration$estimate[1])
    )$estimate$estimate
  )
})
