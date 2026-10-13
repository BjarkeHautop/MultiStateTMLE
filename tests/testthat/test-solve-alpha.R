# solve_alpha() solves the calibration equation
# alpha^theta(P) := (Psi_z(P))^-1(theta) from Rytgaard & van der Laan's
# calibrated stochastic interventions paper (Section 3.3). Lemma 1
# guarantees alpha -> Psi_z^alpha(P) is increasing and concave, an
# assumption the code here relies on.

test_that("solve_alpha() validates its arguments", {
  fun <- function(alpha) list(estimate = c(tmle.est = alpha), eic = alpha)

  expect_error(
    solve_alpha(theta = 1, fun = fun),
    "c_n"
  )

  expect_error(
    solve_alpha(
      theta = 1,
      fun = fun,
      c_n = 1e-4,
      alpha_min = 5,
      alpha_max = 1
    ),
    "0 < alpha_min < alpha_max"
  )

  expect_error(
    solve_alpha(theta = 1, fun = fun, c_n = 1e-4, alpha_init = Inf),
    "alpha_init"
  )
})

test_that("solve_alpha() returns immediately when alpha_init already matches theta", {
  n_calls <- 0
  fun <- function(alpha) {
    n_calls <<- n_calls + 1
    list(estimate = c(tmle.est = alpha, se = 0.1), eic = rep(alpha, 10))
  }

  res <- solve_alpha(
    theta = 1,
    fun = fun,
    c_n = 1e-4,
    alpha_init = 1
  )

  expect_true(res$converged)
  expect_equal(res$alpha.hat, 1)
  expect_equal(n_calls, 1)
})

test_that("solve_alpha() finds alpha for a monotone increasing estimand", {
  # psi(alpha) = alpha, so the root at theta is alpha = theta
  fun <- function(alpha) {
    list(estimate = c(tmle.est = alpha, se = 0.1), eic = rep(alpha - 3, 20))
  }

  res <- solve_alpha(
    theta = 3,
    fun = fun,
    c_n = 1e-4,
    alpha_init = 1
  )

  expect_true(res$converged)
  expect_equal(res$alpha.hat, 3)
  expect_equal(res$eic, rep(res$alpha.hat - 3, 20))
})

test_that("solve_alpha() searches downward when starting above theta", {
  fun <- function(alpha) {
    list(estimate = c(tmle.est = alpha, se = 0.1), eic = rep(0, 5))
  }

  res <- solve_alpha(
    theta = 2,
    fun = fun,
    c_n = 1e-4,
    alpha_init = 50
  )

  expect_true(res$converged)
  expect_equal(res$alpha.hat, 2)
})

test_that("solve_alpha() caches evaluations: every fun call is recorded exactly once", {
  n_calls <- 0
  fun <- function(alpha) {
    n_calls <<- n_calls + 1
    list(estimate = c(tmle.est = alpha, se = 0.1), eic = rep(0, 5))
  }

  res <- solve_alpha(
    theta = 2.5,
    fun = fun,
    c_n = 1e-4,
    alpha_init = 1
  )

  expect_equal(n_calls, nrow(res$grid))
  expect_equal(n_calls, length(unique(res$grid$alpha)))
})

test_that("solve_alpha() solves the calibration equation for a Psi_z-shaped curve", {
  # Psi_z^alpha(P) = L * (1 - exp(-alpha))
  # has analytic inverse
  # alpha = -log(1 - theta / L).
  L <- 0.6
  fun <- function(alpha) {
    psi <- L * (1 - exp(-alpha))
    list(estimate = c(tmle.est = psi, se = 0.05), eic = rep(psi, 10))
  }

  theta <- 0.3
  res <- solve_alpha(
    theta = theta,
    fun = fun,
    c_n = 1e-4,
    alpha_init = 1
  )

  expected_alpha <- -log(1 - theta / L)

  expect_true(res$converged)
  expect_equal(res$alpha.hat, expected_alpha, tolerance = 1e-3)
})

test_that("solve_alpha() reports no bracket found when the estimand is out of range", {
  # psi is constant, so it can never bracket a theta above it
  fun <- function(alpha) {
    list(estimate = c(tmle.est = 1, se = 0.1), eic = rep(0, 5))
  }

  res <- solve_alpha(
    theta = 100,
    fun = fun,
    c_n = 1e-4,
    alpha_init = 1,
    alpha_max = 10
  )

  expect_false(res$converged)
  expect_equal(
    res$alpha.hat,
    res$grid$alpha[which.min(abs(res$grid$psi - 100))]
  )
})
