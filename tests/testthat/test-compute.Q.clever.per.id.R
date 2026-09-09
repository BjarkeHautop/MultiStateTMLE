# Rytgaard, Gerds & van der Laan (2022, Annals of Statistics) define Q as a
# backward recursion in conditional expectations, and the clever covariate
# for a process x as h_t^x = E[Y | jump] - E[Y | no jump] (their Eq. 16-18).

test_that("compute.Q.clever.per.id reproduces the discrete-time cumulative incidence for a single terminal process", {
  # S = 1, single terminal process "D": Q_t = h_t + (1 - h_t) * Q_{t+1}
  # collapses to Q_t = 1 - prod_{s=t}^{T} (1 - h_s).
  h <- c(0.1, 0.2, 0.3)

  dt_id <- data.table::data.table(
    state = c(1L, 1L, 1L),
    state.row.index = c(1L, 1L, 1L),
    P.D.1 = h
  )

  states <- data.table::data.table(state = 1L)

  out <- compute.Q.clever.per.id(
    dt_id = dt_id,
    states = states,
    process.types = list(D = "terminal")
  )

  expected_Q <- rev(1 - cumprod(rev(1 - h)))

  expect_equal(out$Q, expected_Q)
})

test_that("compute.Q.clever.per.id: clever covariates match the analytic recursion", {
  h <- c(0.1, 0.2, 0.3)

  dt_id <- data.table::data.table(
    state = c(1L, 1L, 1L),
    state.row.index = c(1L, 1L, 1L),
    P.D.1 = h
  )
  states <- data.table::data.table(state = 1L)

  out <- compute.Q.clever.per.id(
    dt_id = dt_id,
    states = states,
    process.types = list(D = "terminal")
  )

  Q <- rev(1 - cumprod(rev(1 - h)))

  expect_equal(out$clever.Q.D1, c(1, 1, 1))
  expect_equal(out$clever.Q.D0, c(Q[2], Q[3], Q[3]))
  expect_equal(out$clever.Q.D, out$clever.Q.D1 - out$clever.Q.D0)
})

test_that("compute.Q.clever.per.id: clever covariate for a mediator process matches the jump-vs-no-jump formula", {
  # S = 2, for the state-transition path.
  # Process "Z" occurs and raises the hazard of terminal process "D".
  states <- data.table::data.table(
    state = c(1L, 2L),
    Z = c(0, 1),
    s.Z = c(2L, 2L)
  )

  dt_id <- data.table::data.table(
    state = c(1L, 1L),
    state.row.index = c(1L, 1L),
    P.Z.1 = c(0.2, 0.2),
    P.Z.2 = c(0, 0),
    P.D.1 = c(0.1, 0.1),
    P.D.2 = c(0.3, 0.3)
  )

  out <- compute.Q.clever.per.id(
    dt_id = dt_id,
    states = states,
    process.types = list(Z = "recurrent", D = "terminal"),
    parameter = "D"
  )

  # Q_2(state 1) = 0.1, Q_2(state 2) = 0.3 (terminal hazards).
  # Q_1(state 1) = 0.1 + 0.7 * 0.1 + 0.2 * 0.3 = 0.23
  expect_equal(out$Q, c(0.23, 0.10))

  # h^Z = Q(successor state 2) - Q(stay state 1) = 0.3 - 0.1 = 0.2
  expect_equal(out$clever.Q.Z1, c(0.3, 0.3))
  expect_equal(out$clever.Q.Z0, c(0.1, 0.1))
  expect_equal(out$clever.Q.Z, c(0.2, 0.2))
})

test_that("compute.Q.clever.per.id returns Q = 0 with a warning when no hazard columns are found", {
  dt_id <- data.table::data.table(
    state = c(1L, 1L),
    state.row.index = c(1L, 1L),
    x = c(1, 2)
  )
  states <- data.table::data.table(state = 1L)

  expect_warning(
    out <- compute.Q.clever.per.id(dt_id = dt_id, states = states),
    "No hazard columns found"
  )

  expect_equal(out$Q, c(0, 0))
})

test_that("compute.Q.clever.per.id errors on empty input", {
  states <- data.table::data.table(state = 1L)
  dt_id <- data.table::data.table(
    state = integer(0),
    state.row.index = integer(0),
    P.D.1 = numeric(0)
  )

  expect_error(
    compute.Q.clever.per.id(dt_id = dt_id, states = states),
    "at least one row"
  )
})
