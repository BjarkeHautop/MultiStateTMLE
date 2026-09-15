test_that("tmle.alpha.fun validates its arguments", {
  initial.fit <- make.fixture.initial.fit(n = 30)

  expect_error(
    tmle.alpha.fun(
      initial.fit = initial.fit,
      target = 1,
      tau = 1,
      alpha = 1,
      use.cores = 1
    ),
    "'target'"
  )
  expect_error(
    tmle.alpha.fun(
      initial.fit = initial.fit,
      target = "outcome1",
      tau = 1,
      alpha = -1,
      use.cores = 1
    ),
    "'alpha'"
  )
  expect_error(
    tmle.alpha.fun(
      initial.fit = initial.fit,
      target = "outcome1",
      tau = 1,
      alpha = 1,
      use.cores = 0
    ),
    "'use.cores'"
  )
  expect_error(
    tmle.alpha.fun(
      initial.fit = initial.fit,
      target = "outcome1",
      tau = -1,
      alpha = 1,
      use.cores = 1
    ),
    "'tau'"
  )
  expect_error(
    tmle.alpha.fun(
      initial.fit = initial.fit,
      target = "outcome1",
      tau = 1,
      alpha = 1,
      use.cores = 1,
      verbose = "x"
    ),
    "'verbose'"
  )
})

test_that("tmle.alpha.fun returns the documented estimate structure", {
  initial.fit <- make.fixture.initial.fit(n = 60)

  fit <- tmle.alpha.fun(
    initial.fit = initial.fit,
    target = "outcome1",
    tau = 1,
    alpha = 1,
    use.cores = 1,
    verbose = FALSE,
    output.eic = TRUE
  )

  expect_true(all(c("tmle.est", "se", "g.est") %in% names(fit$estimate)))
  expect_true(is.numeric(fit$eic))
  expect_equal(length(fit$eic), length(unique(initial.fit$tmp.long$id)))
  expect_true(fit$estimate[["tmle.est"]] >= 0 && fit$estimate[["tmle.est"]] <= 1)
})
