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
    tmle_alpha(new_fit, target = "death", alpha = 0.5, z_name = "z")
  )

  # An unknown z.name is ignored, and alpha with it.
  expect_warning(
    old_est <- suppressWarnings(
      tmle.alpha.fun(
        initial.fit = old_fit,
        target = "death",
        alpha = 0.5,
        z.name = "w"
      ),
      classes = "lifecycle_warning_deprecated"
    ),
    "alpha is ignored"
  )
  expect_equal(old_est, tmle_alpha(new_fit, target = "death"))

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

  expect_error(
    estimate.derivative(alpha_hat = 3),
    "calibrate_alpha",
    class = "lifecycle_error_deprecated"
  )
  expect_error(
    estimate.alpha.fun(theta = 2),
    "calibrate_alpha",
    class = "lifecycle_error_deprecated"
  )

  treatment.fit <- fake.tmle.fit(0.7, c(0.1, -0.1))
  placebo.fit <- fake.tmle.fit(0.3, c(0.2, -0.2), a = 0)
  calibrated.fit <- fake.calibrated.fit(0.5, c(0, 0))
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

test_that(".split_control() moves tmle_alpha() settings into control", {
  expect_equal(
    .split_control(list(min_iter = 2, cores = 1)),
    list(cores = 1, control = list(min_iter = 2))
  )
  expect_equal(.split_control(list(cores = 1)), list(cores = 1))
})

test_that("dotted tmle.alpha.fun() arguments are translated for tmle_alpha()", {
  expect_equal(
    .snake_case_args(list(
      min.iter = 2,
      use.cores = 1,
      tau = 3,
      output.weights = 0.5
    )),
    list(min_iter = 2, cores = 1)
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

  expect_equal(
    est,
    tmle_alpha(fixture.fit(30), target = "death", z_name = "z")
  )
})

test_that("deprecated calibration functions translate their targets", {
  withr::local_options(lifecycle_verbosity = "quiet")
  fit <- fixture.fit(100)

  expect_equal(
    calibration.fun(fit, rho = 0.8, target = "death", verbose = FALSE),
    calibrate_alpha(
      fit,
      target = "death",
      z_name = "z",
      value = 0.8,
      scale = "ratio"
    )
  )
  expect_error(
    calibration.fun(fit, rho = 0.8, delta = 0.1, target = "death"),
    "at most one"
  )
  expect_error(
    calibration.fun(fit, rho.1a = 0.8, target = "death"),
    "requires a baseline treatment"
  )

  expect_equal(
    calibration.curve.fun(
      fit,
      target = "death",
      alpha.grid = 1,
      block.size.target = 5,
      verbose = FALSE
    ),
    calibration_curve(
      fit,
      target = "death",
      z_name = "z",
      alpha_grid = 1,
      estimand = "years_lost",
      control = list(min_iter = 1, block_size = 5)
    )
  )
  expect_error(
    calibration.curve.fun(
      fit,
      target = "death",
      block.size.z = 5,
      block.size.target = 10
    ),
    "must be equal"
  )
})
