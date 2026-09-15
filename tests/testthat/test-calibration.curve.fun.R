test_that("calibration.curve.fun validates its arguments", {
  initial.fit <- make.fixture.initial.fit(n = 30)

  expect_error(
    calibration.curve.fun(
      initial.fit = initial.fit,
      alpha.grid = numeric(0),
      tau = 1,
      target = "outcome1",
      z.name = "z",
      use.cores = 1
    ),
    "'alpha.grid'"
  )
  expect_error(
    calibration.curve.fun(
      initial.fit = initial.fit,
      alpha.grid = c(0, 1),
      tau = 1,
      target = 1,
      z.name = "z",
      use.cores = 1
    ),
    "'target'"
  )
  expect_error(
    calibration.curve.fun(
      initial.fit = initial.fit,
      alpha.grid = c(0, 1),
      tau = 1,
      target = "outcome1",
      z.name = "z",
      use.cores = 1,
      min.iter = -1
    ),
    "'min.iter'"
  )
})

test_that("calibration.curve.fun traces the estimand over the alpha grid", {
  initial.fit <- make.fixture.initial.fit(n = 60)

  curve <- calibration.curve.fun(
    initial.fit = initial.fit,
    alpha.grid = c(0, 1),
    tau = 1,
    target = "outcome1",
    z.name = "z",
    use.cores = 1,
    verbose = FALSE
  )

  expect_equal(nrow(curve$estimate), 2)
  expect_equal(curve$estimate$alpha, c(0, 1))
  # alpha = 0 fully suppresses z, so exposure and target estimates should be 0.
  expect_equal(unname(curve$estimate$exposure.est[1]), 0)
})
