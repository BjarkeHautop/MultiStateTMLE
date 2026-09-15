make.fixture.sim.parameters <- function(n = 60, seed = 1405) {
  dt <- make.fixture.dt(n = n, seed = seed)
  prepare.initial(
    dt,
    tau = 1,
    fit.types = make.fixture.fit.types(),
    verbose = FALSE,
    return.parameters.for.simulation = TRUE
  )
}

test_that("sim.from.data validates its arguments", {
  sim.parameters <- make.fixture.sim.parameters()

  expect_error(
    sim.from.data(n = -1, sim.parameters = sim.parameters),
    "'n'"
  )
  expect_error(
    sim.from.data(n = 10, sim.parameters = "x"),
    "'sim.parameters'"
  )
  expect_error(
    sim.from.data(n = 10, sim.parameters = sim.parameters, cens = "x"),
    "'cens'"
  )
  expect_error(
    sim.from.data(n = 10, sim.parameters = sim.parameters, browse = "x"),
    "'browse'"
  )
})

test_that("sim.from.data simulates event-history data with the expected columns", {
  sim.parameters <- make.fixture.sim.parameters()

  set.seed(2)
  new.dt <- sim.from.data(n = 50, sim.parameters = sim.parameters)

  expect_true(data.table::is.data.table(new.dt))
  expect_true(all(c("id", "time", "delta", "L0", "z") %in% names(new.dt)))
  expect_equal(length(unique(new.dt$id)), 50)
  expect_true(all(new.dt$time >= 0))
})

test_that("sim.from.data's alpha.intervention scales the target process's rate", {
  sim.parameters <- make.fixture.sim.parameters(n = 300)

  set.seed(3)
  dt.baseline <- sim.from.data(n = 400, sim.parameters = sim.parameters)
  set.seed(3)
  dt.suppressed <- sim.from.data(
    n = 400,
    sim.parameters = sim.parameters,
    alpha.intervention = list(z = 0.01)
  )

  # Fraction of subjects who ever have z = 1 should drop sharply when the
  # z-intensity is scaled down by a factor of 100.
  frac.z.baseline <- mean(dt.baseline[, max(z), by = "id"][[2]])
  frac.z.suppressed <- mean(dt.suppressed[, max(z), by = "id"][[2]])

  expect_lt(frac.z.suppressed, frac.z.baseline)
})
