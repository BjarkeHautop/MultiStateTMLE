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
