test_that("sim.generic validates its arguments", {
  baseline <- list(L0 = function(N) rnorm(N))
  processes <- list(
    z = list(type = "one.jump", eta = 0.2, nu = 1),
    outcome1 = list(type = "terminal", eta = 0.3, nu = 1),
    censoring = list(type = "censoring", eta = 0.1, nu = 1)
  )
  effects <- list(c("L0", "z", 0.5), c("z", "outcome1", 0.7))

  expect_error(
    sim.generic(baseline = "x", processes = processes, effects = effects, n = 10),
    "'baseline'"
  )
  expect_error(
    sim.generic(baseline = baseline, processes = processes, effects = effects, n = -1),
    "'n'"
  )
  expect_error(
    sim.generic(
      baseline = baseline,
      processes = processes,
      effects = effects,
      n = 10,
      cens = "x"
    ),
    "'cens'"
  )
  expect_error(
    sim.generic(
      baseline = baseline,
      processes = processes,
      effects = effects,
      n = 10,
      alpha.intervention = "x"
    ),
    "'alpha.intervention'"
  )
})

test_that("sim.generic produces reproducible event-history data with the expected columns", {
  dt1 <- make.fixture.dt(n = 30, seed = 7)
  dt2 <- make.fixture.dt(n = 30, seed = 7)

  expect_true(data.table::is.data.table(dt1))
  expect_true(all(c("id", "time", "delta", "L0", "z") %in% names(dt1)))
  expect_true(nrow(dt1) > 0)
  expect_equal(length(unique(dt1$id)), 30)
  expect_equal(dt1, dt2)
})

test_that("sim.generic's baseline.intervention fixes the covariate for everyone", {
  baseline <- list(L0 = function(N) rnorm(N))
  processes <- list(
    z = list(type = "one.jump", eta = 0.2, nu = 1),
    outcome1 = list(type = "terminal", eta = 0.3, nu = 1),
    censoring = list(type = "censoring", eta = 0.1, nu = 1)
  )
  effects <- list(c("L0", "z", 0.5), c("z", "outcome1", 0.7))

  set.seed(1)
  dt <- sim.generic(
    baseline = baseline,
    processes = processes,
    effects = effects,
    n = 30,
    baseline.intervention = list(L0 = 0)
  )

  expect_true(all(dt$L0 == 0))
})

test_that("sim.generic's cens = 0 disables the censoring process", {
  baseline <- list(L0 = function(N) rnorm(N))
  processes <- list(
    z = list(type = "one.jump", eta = 0.2, nu = 1),
    outcome1 = list(type = "terminal", eta = 0.3, nu = 1),
    censoring = list(type = "censoring", eta = 0.1, nu = 1)
  )
  effects <- list(c("L0", "z", 0.5), c("z", "outcome1", 0.7))

  set.seed(1)
  dt <- sim.generic(
    baseline = baseline,
    processes = processes,
    effects = effects,
    n = 100,
    cens = 0
  )

  expect_equal(sum(dt$delta == 0), 0)
})
