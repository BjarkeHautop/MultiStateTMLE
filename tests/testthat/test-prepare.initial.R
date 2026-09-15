test_that("prepare.initial validates its arguments", {
  dt <- make.fixture.dt(n = 20)

  expect_error(
    prepare.initial(dt = "x", tau = 1),
    "'dt'"
  )
  expect_error(
    prepare.initial(dt = dt, tau = -1),
    "'tau'"
  )
  expect_error(
    prepare.initial(dt = dt[, .(id, time)], tau = 1),
    "subset"
  )
  expect_error(
    prepare.initial(dt = dt, tau = 1, max.count = -1),
    "'max.count'"
  )
  expect_error(
    prepare.initial(dt = dt, tau = 1, verbose = "x"),
    "'verbose'"
  )
})

test_that("prepare.initial returns the documented fit object structure", {
  initial.fit <- make.fixture.initial.fit(n = 40)

  expect_true(is.list(initial.fit))
  expect_true(all(
    c(
      "tmp.long",
      "depend.matrix",
      "process.names",
      "process.types",
      "process.deltas",
      "cens.process.id",
      "at.risks"
    ) %in%
      names(initial.fit)
  ))

  expect_true(data.table::is.data.table(initial.fit$tmp.long))
  expect_equal(initial.fit$process.names, c("z", "outcome1", "censoring"))
  expect_equal(
    length(unique(initial.fit$tmp.long$id)),
    40
  )
})

test_that("prepare.initial with return.parameters.for.simulation = TRUE returns simulation parameters", {
  dt <- make.fixture.dt(n = 40)
  params <- prepare.initial(
    dt,
    tau = 1,
    fit.types = make.fixture.fit.types(),
    verbose = FALSE,
    return.parameters.for.simulation = TRUE
  )

  expect_true(all(c("z", "outcome1", "censoring", "model.structure") %in% names(params)))
  expect_true("weibull.parameters" %in% names(params$z))
  expect_equal(params$n.subjects, 40)
})

test_that("prepare.initial summarizes a categorical baseline variable", {
  dt <- data.table::copy(make.fixture.dt(n = 40))
  dt[, G := sample(c("a", "b"), 1), by = "id"]

  fit.types <- make.fixture.fit.types()
  fit.types$z$model <- "Surv(tstart, tstop, delta == 2)~L0+G"

  params <- prepare.initial(
    dt,
    tau = 1,
    fit.types = fit.types,
    verbose = FALSE,
    return.parameters.for.simulation = TRUE
  )

  expect_equal(params$baseline.summary$G$type, "character")
  expect_setequal(names(params$baseline.summary$G$frequencies), c("a", "b"))
})

test_that("prepare.initial converts a non-numeric time column with a warning", {
  dt <- make.fixture.dt(n = 20)
  dt[, time := as.character(time)]

  expect_warning(
    initial.fit <- prepare.initial(
      dt,
      tau = 1,
      fit.types = make.fixture.fit.types(),
      verbose = FALSE
    ),
    "not numeric"
  )
  expect_true(is.numeric(initial.fit$tmp.long$time))
})

test_that("prepare.initial fits a treatment model and computes propensities under intervention 'a'", {
  dt <- data.table::copy(make.fixture.dt(n = 40))
  dt[, A0 := rbinom(1, 1, 0.5), by = "id"]

  initial.fit <- prepare.initial(
    dt,
    tau = 1,
    fit.types = make.fixture.fit.types(),
    fit.treatment = list(model = "A0~L0", fit = "glm"),
    a = c(0, 1),
    verbose = FALSE
  )

  expect_true("pi.A0.1" %in% names(initial.fit$tmp.long))
  expect_true(all(
    initial.fit$tmp.long$pi.A0.1 >= 0 & initial.fit$tmp.long$pi.A0.1 <= 1
  ))
})

test_that("prepare.initial reports unsupported treatment fit methods", {
  dt <- data.table::copy(make.fixture.dt(n = 20))
  dt[, A0 := rbinom(1, 1, 0.5), by = "id"]

  expect_output(
    prepare.initial(
      dt,
      tau = 1,
      fit.types = make.fixture.fit.types(),
      fit.treatment = list(model = "A0~L0", fit = "other"),
      verbose = FALSE
    ),
    "need to incorporate"
  )
})

test_that("prepare.initial prints fit diagnostics when verbose=TRUE", {
  dt <- make.fixture.dt(n = 20)

  expect_output(
    prepare.initial(
      dt,
      tau = 1,
      fit.types = make.fixture.fit.types(),
      verbose = TRUE
    )
  )
})

test_that("prepare.initial's HAL setup identifies the at-risk process for its own delta", {
  # Regression test for a `name.hal` indexing bug that used to select every
  # process name instead of just the HAL one and error before reaching the
  # HAL fit itself.
  dt <- make.fixture.dt(n = 30)
  fit.types <- make.fixture.fit.types()
  fit.types$z$fit <- "hal"

  expect_error(
    prepare.initial(dt, tau = 1, fit.types = fit.types, verbose = FALSE),
    "could not find function \"fit.hal\"",
    fixed = TRUE
  )
})

test_that("HAL fitting beyond fit.hal() is untested", {
  skip(
    "fit.hal() is called by prepare.initial() but is not defined anywhere."
  )
})
