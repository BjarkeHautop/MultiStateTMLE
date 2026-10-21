test_that("deprecated estimation functions warn and match their replacements", {
  withr::local_options(lifecycle_verbosity = "warning")

  old_fit_types <- list(
    z = list(
      model = "Surv(tstart, tstop, delta == 2)~L0",
      fit = "cox",
      at.risk = function(dt) (dt[["z"]] == 0)
    ),
    death = list(
      model = "Surv(tstart, tstop, delta == 1)~L0+z",
      fit = "cox"
    ),
    censoring = list(model = "Surv(tstart, tstop, delta == 0)~L0", fit = "cox")
  )
  expect_warning(
    old_fit <- prepare.initial(
      illness_death,
      tau = 1,
      fit.types = old_fit_types,
      use.cores = 1
    ),
    "fit_initial",
    class = "lifecycle_warning_deprecated"
  )
  new_fit <- fit_initial(
    illness_death,
    fit_types = illness.death.fit.types,
    tau = 1
  )
  expect_equal(old_fit$tmp_long, new_fit$tmp_long)

  expect_warning(
    old_est <- tmle.alpha.fun(
      initial.fit = old_fit,
      target = "death",
      alpha = 0.5,
      use.cores = 1
    ),
    "tmle_alpha",
    class = "lifecycle_warning_deprecated"
  )
  expect_equal(
    old_est,
    tmle_alpha(new_fit, target = "death", alpha = 0.5)
  )

  expect_error(
    .prepare_initial(
      illness_death,
      tau = 1,
      fit.types = old_fit_types,
      return.parameters.for.simulation = TRUE
    ),
    "sim_model_from_fits",
    class = "lifecycle_error_deprecated"
  )
})

test_that("deprecated helper functions warn and match their replacements", {
  withr::local_options(lifecycle_verbosity = "warning")

  fun <- function(alpha, parameter) list(estimate = c(tmle_est = alpha^2))
  expect_warning(
    d <- estimate.derivative(alpha_hat = 3, fun = fun, h = 0.01),
    "estimate_derivative",
    class = "lifecycle_warning_deprecated"
  )
  expect_equal(d, 6)

  fun <- function(alpha) list(estimate = c(tmle_est = alpha, se = 0.1), eic = 0)
  expect_warning(
    res <- estimate.alpha.fun(theta = 2, fun = fun, c_n = 1e-4),
    "solve_alpha",
    class = "lifecycle_warning_deprecated"
  )
  expect_equal(res, solve_alpha(fun, theta = 2, c_n = 1e-4))

  treatment.fit <- list(estimate = c(tmle_est = 0.7), eic = c(0.1, -0.1))
  placebo.fit <- list(estimate = c(tmle_est = 0.3), eic = c(0.2, -0.2))
  calibrated.fit <- list(estimate = c(target_est = 0.5), eic = c(0, 0))
  expect_warning(
    contrasts <- make.calibrated.contrasts(
      treatment.fit = treatment.fit,
      placebo.fit = placebo.fit,
      calibrated.fit = calibrated.fit
    ),
    "calibrated_contrasts",
    class = "lifecycle_warning_deprecated"
  )
  expect_equal(
    contrasts,
    calibrated_contrasts(treatment.fit, placebo.fit, calibrated.fit)
  )

  dt_id <- data.table::data.table(
    state = 1L,
    state.row.index = 1L,
    P.D.1 = 0.1
  )
  states <- data.table::data.table(state = 1L)
  expect_warning(
    compute.Q.clever.per.id(
      dt_id,
      states,
      process.types = list(D = "terminal")
    ),
    class = "lifecycle_warning_deprecated"
  )
})

test_that("dotted tmle.alpha.fun() arguments are translated for tmle_alpha()", {
  expect_equal(
    .snake_case_args(list(
      min.iter = 2,
      use.cores = 1,
      tau = 3,
      output.weights = 0.5
    )),
    list(min_iter = 2, cores = 1, output_weights = 0.5)
  )
  expect_equal(.snake_case_args(list()), list())
})

test_that("deprecated tmle.alpha.fun() builds the initial fit from dt", {
  withr::local_options(lifecycle_verbosity = "quiet")

  est <- tmle.alpha.fun(
    dt = fixture.data(30),
    fit.types = list(
      z = list(
        model = "Surv(tstart, tstop, delta == 2)~L0",
        at.risk = function(dt) (dt[["z"]] == 0)
      ),
      death = "Surv(tstart, tstop, delta == 1)~L0+z",
      censoring = "Surv(tstart, tstop, delta == 0)~L0"
    ),
    target = "death",
    tau = 1,
    use.cores = 1
  )

  expect_equal(est, tmle_alpha(fixture.fit(30), target = "death"))
})
