test_that("print methods summarise the result objects", {
  fit <- fixture.fit(30)
  expect_output(print(fit), "Initial fit of 30 individuals up to tau = 1")
  expect_output(print(fit), "one_jump")

  est <- tmle_alpha(fit, target = "death", alpha = 0.5, z_name = "z")
  expect_output(print(est), "TMLE estimate under alpha = 0.5 on z")
  expect_output(print(est), "converged after")
  expect_output(
    print(tmle_alpha(fit, target = "death", one_step = TRUE)),
    "One-step estimate without intervention"
  )

  curve <- calibration_curve(
    fit,
    target = "death",
    z_name = "z",
    alpha_grid = 1
  )
  expect_output(print(curve), "Calibration curve over 1 interventions on z")
})

test_that("print shows the treatment model of an initial fit", {
  fit <- fixture.fit(30)
  fit$treatment_model <- A0 ~ L0
  fit$a <- 0:1

  expect_output(print(fit), "Treatment model: A0 ~ L0, prepared for a = 0, 1")
})

test_that("coef() and confint() extract the estimates by target", {
  est <- fake.tmle.fit(0.3, 0)
  est$estimate$se <- 0.1

  expect_equal(coef(est), c(death = 0.3))
  expect_equal(
    confint(est),
    matrix(
      0.3 + c(-1, 1) * qnorm(0.975) * 0.1,
      nrow = 1,
      dimnames = list("death", c("2.5 %", "97.5 %"))
    )
  )
  expect_equal(
    colnames(confint(est, level = 0.9)),
    c("5 %", "95 %")
  )

  cal <- fake.calibrated.fit(c(0.2, 0.4), list(0, 0), target = c("z", "death"))
  cal$estimate$se <- c(0.1, 0.2)

  expect_equal(coef(cal), c(z = 0.2, death = 0.4))
  expect_equal(rownames(confint(cal, "death")), "death")
  expect_equal(
    confint(cal, "death")[1, ],
    0.4 + c(-1, 1) * qnorm(0.975) * 0.2,
    ignore_attr = TRUE
  )
  expect_error(confint(cal, "stroke"), "stroke")
})
