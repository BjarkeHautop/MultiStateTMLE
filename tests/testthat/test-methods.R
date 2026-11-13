test_that("print methods summarise the result objects", {
  fit <- fixture.fit(30)
  expect_output(print(fit), "Initial fit of 30 individuals up to tau = 1")
  expect_output(print(fit), "one_jump")

  est <- tmle_estimate(fit, risk(death), scale_intensity(illness = 0.5))
  expect_output(print(est), "TMLE estimate of risk(death)", fixed = TRUE)
  expect_output(
    print(est),
    "under scale_intensity(illness = 0.5)",
    fixed = TRUE
  )
  expect_output(
    print(est),
    .format_targeting(est$checks, "death"),
    fixed = TRUE
  )
  expect_output(
    print(
      tmle_estimate(fit, risk(death), control = list(one_step = TRUE))
    ),
    "One-step estimate of risk(death)\nwithout intervention",
    fixed = TRUE
  )

  curve <- intervention_curve(fit, risk(death), risk(illness), alpha = 1)
  expect_output(
    print(curve),
    "over 1 scalings of the intensity of illness"
  )
})

test_that("print shows the treatment model of an initial fit", {
  fit <- fixture.fit(30)
  fit$treatment_model <- A0 ~ L0
  fit$arms <- 0:1

  expect_output(
    print(fit),
    "Treatment model: A0 ~ L0, prepared for arms a = 0, 1"
  )
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

  cal <- fake.calibrated.fit(
    c(0.2, 0.4),
    list(0, 0),
    target = c("illness", "death")
  )
  cal$estimate$se <- c(0.1, 0.2)

  expect_equal(coef(cal), c(illness = 0.2, death = 0.4))
  expect_equal(rownames(confint(cal, "death")), "death")
  expect_equal(
    confint(cal, "death")[1, ],
    0.4 + c(-1, 1) * qnorm(0.975) * 0.2,
    ignore_attr = TRUE
  )
  expect_error(confint(cal, "stroke"), "stroke")
})

test_that("the targeting summary counts updates, not iterations", {
  checks <- function(iter, converged) {
    list(iter = c(death = iter), converged = c(death = converged))
  }
  expect_equal(
    .format_targeting(checks(1, TRUE), "death"),
    "no targeting needed, as it already solves the EIC equation"
  )
  expect_equal(
    .format_targeting(checks(2, TRUE), "death"),
    "converged after 1 targeting step"
  )
  expect_equal(
    .format_targeting(checks(4, TRUE), "death"),
    "converged after 3 targeting steps"
  )
  expect_equal(
    .format_targeting(checks(10, FALSE), "death"),
    "did not converge after 10 targeting steps"
  )
})
