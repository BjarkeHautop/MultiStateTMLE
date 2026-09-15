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

test_that("tmle.alpha.fun builds initial.fit internally from 'dt' when not supplied", {
  dt <- make.fixture.dt(n = 30)

  fit <- tmle.alpha.fun(
    dt = dt,
    fit.types = make.fixture.fit.types(),
    target = "outcome1",
    tau = 1,
    alpha = 1,
    z.name = "z",
    use.cores = 1,
    verbose = FALSE
  )

  expect_true(fit$estimate[["tmle.est"]] >= 0 && fit$estimate[["tmle.est"]] <= 1)
})

test_that("tmle.alpha.fun warns and ignores alpha when z.name matches no process", {
  initial.fit <- make.fixture.initial.fit(n = 30)

  expect_warning(
    fit <- tmle.alpha.fun(
      initial.fit = initial.fit,
      target = "outcome1",
      tau = 1,
      alpha = 2,
      z.name = "not.a.process",
      use.cores = 1,
      verbose = FALSE
    ),
    "No 'z' process found"
  )
  expect_true(fit$estimate[["tmle.est"]] >= 0 && fit$estimate[["tmle.est"]] <= 1)
})

test_that("tmle.alpha.fun supports only.first, target.only.in.state, target.by.state, and one.step", {
  initial.fit <- make.fixture.initial.fit(n = 30)

  fit.only.first <- tmle.alpha.fun(
    initial.fit = initial.fit,
    target = "z",
    tau = 1,
    alpha = 1,
    only.first = TRUE,
    use.cores = 1,
    verbose = FALSE
  )
  expect_true(is.numeric(fit.only.first$estimate[["tmle.est"]]))

  fit.in.state <- tmle.alpha.fun(
    initial.fit = initial.fit,
    target = "outcome1",
    tau = 1,
    alpha = 1,
    target.only.in.state = function(states) states[["z"]] == 0,
    use.cores = 1,
    verbose = FALSE
  )
  expect_true(is.numeric(fit.in.state$estimate[["tmle.est"]]))

  fit.by.state <- tmle.alpha.fun(
    initial.fit = initial.fit,
    target = "outcome1",
    tau = 1,
    alpha = 1,
    target.by.state = TRUE,
    use.cores = 1,
    verbose = FALSE
  )
  expect_true(is.numeric(fit.by.state$estimate[["tmle.est"]]))

  fit.one.step <- tmle.alpha.fun(
    initial.fit = initial.fit,
    target = "outcome1",
    tau = 1,
    alpha = 1,
    one.step = TRUE,
    use.cores = 1,
    verbose = FALSE
  )
  expect_true("one.step.est" %in% names(fit.one.step$estimate))
})

test_that("tmle.alpha.fun truncates clever weights and reports weight/convergence summaries", {
  initial.fit <- make.fixture.initial.fit(n = 30)

  fit <- tmle.alpha.fun(
    initial.fit = initial.fit,
    target = "outcome1",
    tau = 1,
    alpha = 3,
    truncate.weights = 1.5,
    output.weights = c(0.5),
    output.convergence = TRUE,
    use.cores = 1,
    verbose = FALSE
  )

  expect_lte(unname(fit$weights[["q50"]]), 1.5 + 1e-8)
  expect_true("no.truncated" %in% names(fit$weights))
  expect_true("cens.weights" %in% names(fit))
  expect_true(all(c("iter", "eic.solved.at", "converged") %in% names(fit$convergence)))
})

test_that("tmle.alpha.fun computes weights and estimates under a baseline intervention 'a'", {
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

  fit <- tmle.alpha.fun(
    initial.fit = initial.fit,
    target = "outcome1",
    tau = 1,
    alpha = 1,
    a = 1,
    output.a.weights = c(0.5),
    use.cores = 1,
    verbose = FALSE
  )
  expect_true(is.numeric(fit$a.weights))
  expect_true(fit$estimate[["tmle.est"]] >= 0 && fit$estimate[["tmle.est"]] <= 1)

  expect_error(
    tmle.alpha.fun(
      initial.fit = initial.fit,
      target = "outcome1",
      tau = 1,
      alpha = 1,
      a = 5,
      use.cores = 1,
      verbose = FALSE
    ),
    "'a' must match the intervention"
  )
})

test_that("tmle.alpha.fun's alpha.list option reproduces the scalar alpha result for a constant intervention", {
  # alpha.list is the (still-experimental, per its "not made yet" comment in
  # R/tmle.alpha.fun.R) way to pass a per-time/covariate intervention
  # function instead of a single scalar alpha.
  initial.fit <- make.fixture.initial.fit(n = 30)

  fit.scalar <- tmle.alpha.fun(
    initial.fit = initial.fit,
    target = "outcome1",
    tau = 1,
    alpha = 2,
    z.name = "z",
    use.cores = 1,
    verbose = FALSE
  )

  expect_message(
    fit.list <- tmle.alpha.fun(
      initial.fit = initial.fit,
      target = "outcome1",
      tau = 1,
      alpha.list = list(z = function(time) rep(2, length(time))),
      use.cores = 1,
      verbose = FALSE
    ),
    "using new alpha option"
  )

  expect_equal(
    fit.list$estimate[["tmle.est"]],
    fit.scalar$estimate[["tmle.est"]]
  )

  fit.list.by.state <- tmle.alpha.fun(
    initial.fit = initial.fit,
    target = "outcome1",
    tau = 1,
    alpha.list = list(z = function(time) rep(2, length(time))),
    target.by.state = TRUE,
    use.cores = 1,
    verbose = FALSE
  )
  expect_equal(
    fit.list.by.state$estimate[["tmle.est"]],
    fit.scalar$estimate[["tmle.est"]]
  )
})
