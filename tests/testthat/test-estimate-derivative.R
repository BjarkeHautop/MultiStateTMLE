test_that("estimate_derivative() computes a symmetric finite difference", {
  fun <- function(alpha) fake.estimate(alpha^2)

  expect_equal(estimate_derivative(fun, alpha = 3, h = 0.01), 6)
})

test_that("estimate_derivative() shrinks h when it exceeds alpha", {
  calls <- numeric()
  fun <- function(alpha) {
    calls <<- c(calls, alpha)
    fake.estimate(alpha^2)
  }

  estimate_derivative(fun, alpha = 1, h = 10)

  expect_equal(calls, c(1 + 3 / 4, 1 - 3 / 4))
})

test_that("estimate_derivative() passes extra arguments through to fun", {
  fun <- function(alpha, offset) fake.estimate(alpha + offset)

  expect_equal(estimate_derivative(fun, alpha = 5, h = 0.5, offset = 100), 1)
})
