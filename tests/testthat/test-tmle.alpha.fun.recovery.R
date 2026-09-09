test_that("tmle.alpha.fun recovers the true risk with no intervention (alpha=1)", {
  fit <- simulate.recovery.fit(alpha = 1, n = 800, seed = 42)
  truth <- analytic.truth(alpha = 1)

  expect_lt(
    abs(fit$estimate[["tmle.est"]] - truth) / fit$estimate[["se"]],
    4
  )
})

test_that("tmle.alpha.fun recovers the true risk under a shape intervention (alpha=0.5)", {
  fit <- simulate.recovery.fit(alpha = 0.5, n = 800, seed = 42)
  truth <- analytic.truth(alpha = 0.5)

  expect_lt(
    abs(fit$estimate[["tmle.est"]] - truth) / fit$estimate[["se"]],
    4
  )
})
