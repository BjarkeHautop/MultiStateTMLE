test_that("calibration.fun validates its arguments", {
  initial.fit <- make.fixture.initial.fit(n = 30)

  expect_error(
    calibration.fun(
      initial.fit = initial.fit,
      theta = 0.2,
      rho = 0.5,
      delta = 0.1,
      tau = 1,
      target = "outcome1",
      z.name = "z",
      use.cores = 1
    ),
    "Only one of 'rho', 'delta', or 'rho.1a'"
  )
  expect_error(
    calibration.fun(
      initial.fit = initial.fit,
      theta = "x",
      tau = 1,
      target = "outcome1",
      z.name = "z",
      use.cores = 1
    ),
    "'theta'"
  )
  expect_error(
    calibration.fun(
      initial.fit = initial.fit,
      theta = 0.2,
      tau = 1,
      target = c("outcome1", "z"),
      z.name = "z",
      use.cores = 1
    ),
    "'target'"
  )
  expect_error(
    calibration.fun(
      initial.fit = initial.fit,
      theta = 0.2,
      tau = 1,
      target = "outcome1",
      z.name = "z",
      use.cores = 0
    ),
    "'use.cores'"
  )
})

test_that("calibration.fun finds an alpha calibrating the z-process to theta", {
  initial.fit <- make.fixture.initial.fit(n = 400)

  cal <- calibration.fun(
    initial.fit = initial.fit,
    theta = 0.2,
    tau = 1,
    target = "outcome1",
    z.name = "z",
    use.cores = 1,
    verbose = FALSE
  )

  expect_true(all(c("alpha.est", "target.est", "theta") %in% names(cal$estimate)))
  expect_equal(unname(cal$estimate[["theta"]]), 0.2)
  expect_true(cal$checks[["est.alpha.converged"]] == 1)

  # The calibrated alpha should make the z-process's own estimand match theta.
  z.at.alpha.hat <- tmle.alpha.fun(
    initial.fit = initial.fit,
    target = "z",
    tau = 1,
    alpha = cal$estimate[["alpha.est"]],
    z.name = "z",
    use.cores = 1,
    verbose = FALSE
  )
  expect_lt(
    abs(z.at.alpha.hat$estimate[["tmle.est"]] - cal$estimate[["theta"]]),
    0.02
  )
})
