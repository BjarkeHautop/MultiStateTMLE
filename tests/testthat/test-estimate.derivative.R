test_that("estimate.derivative computes a symmetric finite difference", {
  fun <- function(alpha, parameter, ...) {
    list(estimate = c(tmle.est = alpha^2))
  }

  d <- estimate.derivative(alpha_hat = 3, fun = fun, parameter = "z", h = 0.01)

  expect_equal(unname(d), 6)
})

test_that("estimate.derivative shrinks h when it exceeds alpha_hat", {
  calls <- list()
  fun <- function(alpha, parameter, ...) {
    calls[[length(calls) + 1]] <<- alpha
    list(estimate = c(tmle.est = alpha^2))
  }

  alpha_hat <- 1
  h <- 10
  estimate.derivative(alpha_hat = alpha_hat, fun = fun, parameter = "z", h = h)

  used_h <- (calls[[1]] - alpha_hat)
  expect_equal(used_h, alpha_hat * 3 / 4)
})

test_that("estimate.derivative picks up one.step.est when tmle.est is absent", {
  fun <- function(alpha, parameter, ...) {
    list(estimate = c(one.step.est = alpha^2))
  }

  d <- estimate.derivative(alpha_hat = 2, fun = fun, parameter = "z", h = 0.01)

  expect_equal(unname(d), 4)
})

test_that("estimate.derivative passes extra arguments through to fun", {
  fun <- function(alpha, parameter, offset) {
    list(estimate = c(tmle.est = alpha + offset))
  }

  d <- estimate.derivative(
    alpha_hat = 5,
    fun = fun,
    parameter = "z",
    h = 0.5,
    offset = 100
  )

  expect_equal(unname(d), 1)
})
