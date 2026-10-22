test_that("tmle_alpha() returns the estimate and diagnostics", {
  fit <- fit_initial(
    illness_death,
    fit_types = illness.death.fit.types,
    tau = 1
  )

  est <- tmle_alpha(fit, target = "death", alpha = 0.5, z_name = "z")

  expect_s3_class(est, "tmle_alpha")
  expect_named(est$estimate, c("target", "estimate", "se", "lower", "upper"))
  expect_equal(est$estimate$target, "death")
  expect_true(est$estimate$estimate > 0 && est$estimate$estimate < 1)
  expect_equal(
    c(est$estimate$lower, est$estimate$upper),
    est$estimate$estimate + c(-1, 1) * qnorm(0.975) * est$estimate$se
  )
  expect_equal(est$method, "tmle")
  expect_equal(est$alpha, 0.5)
  expect_equal(est$z_name, "z")
  expect_named(est$eic, "death")
  expect_length(est$eic$death, length(unique(illness_death$id)))
  expect_named(
    est$checks,
    c("iter", "eic_solved_at", "converged", "weights", "n_truncated")
  )
  expect_true(est$checks$converged)
  weights <- est$checks$weights
  expect_equal(
    dimnames(weights),
    list(
      c("intervention", "censoring", "treatment", "total"),
      c("q50", "q90", "q99", "q100")
    )
  )
  expect_true(all(weights >= 0))
  # Without a baseline intervention, the treatment weights are 1.
  expect_equal(unname(weights["treatment", ]), rep(1, 4))
  expect_equal(est$checks$n_truncated, 0)
  # The maximum of a product is at most the product of the maxima.
  expect_lte(
    weights["total", "q100"],
    prod(weights[c("intervention", "censoring"), "q100"]) + 1e-8
  )

  one_step <- tmle_alpha(
    fit,
    target = "death",
    alpha = 0.5,
    z_name = "z",
    control = list(one_step = TRUE)
  )
  expect_equal(one_step$method, "one_step")
  expect_named(one_step$estimate, names(est$estimate))
  expect_equal(one_step$initial_estimate, est$initial_estimate)
})

test_that("tmle_alpha() matches process names exactly", {
  fit <- fixture.fit(30)

  expect_error(tmle_alpha(fit, target = "DEATH"), "target must be one of")
  expect_error(
    tmle_alpha(fit, target = "death", z_name = "Z", alpha = 0.5),
    "z_name must be one of"
  )
})

test_that("tmle_alpha() validates its arguments", {
  fit <- fit_initial(
    illness_death,
    fit_types = illness.death.fit.types,
    tau = 1
  )

  expect_error(tmle_alpha(list(), target = "death"), "initial_fit")
  expect_error(tmle_alpha(fit, target = "stroke"), "target must be one of")
  expect_error(tmle_alpha(fit, target = "death", alpha = -1), "alpha")
  expect_error(
    tmle_alpha(fit, target = "death", a = 1),
    "a requires a treatment_model in fit_initial()",
    fixed = TRUE
  )
  expect_error(
    tmle_alpha(fit, target = "death", alpha = 0.5),
    "z_name is required"
  )
  expect_error(
    tmle_alpha(fit, target = "death", alpha = 0.5, z_name = "w"),
    "z_name must be one of"
  )
  # z_name is optional without an intervention, but still validated.
  expect_equal(
    tmle_alpha(fit, target = "death")$estimate,
    tmle_alpha(fit, target = "death", z_name = "z")$estimate
  )
  expect_error(
    tmle_alpha(fit, target = "death", z_name = "w"),
    "z_name must be one of"
  )
  expect_error(
    tmle_alpha(
      fit,
      target = "death",
      alpha = 2,
      alpha_list = list(z = function(time) rep(2, length(time)))
    ),
    "either alpha_list"
  )
  expect_error(
    tmle_alpha(fit, target = "death", control = list(foo = 1)),
    "control"
  )
  expect_error(
    tmle_alpha(fit, target = "death", control = list(conv_const = 0)),
    "control\\$conv_const"
  )
  expect_error(
    tmle_alpha(fit, target = "death", control = list(truncate_weights = 0)),
    "control\\$truncate_weights"
  )
})

test_that("tmle_alpha() supports only_first, target_only_in_state, target_by_state", {
  fit <- fixture.fit(30)

  expect_true(is.numeric(
    tmle_alpha(
      fit,
      target = "z",
      z_name = "z",
      only_first = TRUE
    )$estimate$estimate
  ))
  expect_true(is.numeric(
    tmle_alpha(
      fit,
      target = "death",
      target_only_in_state = function(states) states[["z"]] == 0
    )$estimate$estimate
  ))
  expect_true(is.numeric(
    tmle_alpha(
      fit,
      target = "death",
      control = list(target_by_state = TRUE)
    )$estimate$estimate
  ))
})

test_that("tmle_alpha() truncates clever weights", {
  est <- tmle_alpha(
    fixture.fit(30),
    target = "death",
    alpha = 3,
    z_name = "z",
    control = list(truncate_weights = 1.5)
  )

  expect_lte(est$checks$weights["intervention", "q100"], 1.5 + 1e-8)
  expect_gt(est$checks$n_truncated, 0)
})

test_that("tmle_alpha() estimates under a baseline intervention 'a'", {
  fit <- fit_initial(
    fixture.data.treated(),
    fit_types = illness.death.fit.types,
    tau = 1,
    treatment_model = A0 ~ L0,
    arms = 0:1
  )

  est <- tmle_alpha(fit, target = "death", a = 1)
  # The treatment weights are 1 / P(A0 = 1 | L0) > 1 for the treated, and
  # separate from the censoring weights.
  expect_gt(est$checks$weights["treatment", "q100"], 1)
  expect_equal(
    est$checks$weights["censoring", ],
    tmle_alpha(fit, target = "death", a = 0)$checks$weights["censoring", ]
  )
  expect_true(
    est$estimate$estimate >= 0 && est$estimate$estimate <= 1
  )

  expect_error(tmle_alpha(fit, target = "death", a = 5), "a")
})

test_that("tmle_alpha()'s alpha_list reproduces a constant scalar alpha", {
  fit <- fixture.fit(30)
  constant <- list(z = function(time) rep(2, length(time)))

  scalar <- tmle_alpha(fit, target = "death", alpha = 2, z_name = "z")
  expect_equal(
    tmle_alpha(fit, target = "death", alpha_list = constant)$estimate,
    scalar$estimate
  )
  expect_equal(
    tmle_alpha(
      fit,
      target = "death",
      alpha_list = constant,
      control = list(target_by_state = TRUE)
    )$estimate$estimate,
    scalar$estimate$estimate
  )
})

test_that("tmle_alpha() estimates the expected years lost", {
  fit <- fixture.fit(30)

  risk <- tmle_alpha(fit, target = "death")
  years_lost <- tmle_alpha(fit, target = "death", estimand = "years_lost")

  expect_equal(years_lost$estimand, "years_lost")
  # The integral of the risk up to tau = 1 is below the risk at tau.
  expect_lt(years_lost$estimate$estimate, risk$estimate$estimate)
  expect_gt(years_lost$estimate$estimate, 0)
  expect_output(print(years_lost), "estimate of years lost")
  expect_error(
    tmle_alpha(fit, target = "death", estimand = "rate"),
    "should be one of"
  )
  expect_error(
    tmle_alpha(fit, target = "death", control = list(block_size = 0)),
    "control\\$block_size"
  )
})

test_that("tmle_alpha() only accepts a values prepared by fit_initial()", {
  fit <- fit_initial(
    fixture.data.treated(),
    fit_types = illness.death.fit.types,
    tau = 1,
    treatment_model = A0 ~ L0,
    arms = 1
  )

  expect_equal(tmle_alpha(fit, target = "death", a = 1)$a, 1)
  expect_error(
    tmle_alpha(fit, target = "death", a = 0),
    "a must be one of the values prepared by fit_initial(): 1.",
    fixed = TRUE
  )
  expect_error(
    calibrate_alpha(
      fit,
      target = "death",
      z_name = "z",
      value = 1,
      scale = "ratio",
      a = 1,
      reference_a = 0
    ),
    "reference_a must be one of"
  )
})

test_that("tmle_alpha() estimates several targets separately", {
  fit <- fixture.fit(30)

  both <- tmle_alpha(fit, target = c("death", "z"), z_name = "z", alpha = 2)
  death <- tmle_alpha(fit, target = "death", z_name = "z", alpha = 2)
  z <- tmle_alpha(fit, target = "z", z_name = "z", alpha = 2)

  expect_equal(both$estimate, rbind(death$estimate, z$estimate))
  expect_equal(
    both$initial_estimate,
    c(death$initial_estimate, z$initial_estimate)
  )
  expect_equal(both$eic, c(death$eic, z$eic))
  expect_named(both$checks$converged, c("death", "z"))
  expect_equal(both$checks$weights, death$checks$weights)
  expect_equal(coef(both), c(death = coef(death), z = coef(z)), ignore_attr = TRUE)
  expect_output(print(both), "z: Initial estimate")
  expect_error(tmle_alpha(fit, target = c("death", "death")), "target")
})
