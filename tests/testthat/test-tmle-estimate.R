test_that("tmle_estimate() returns the estimate and diagnostics", {
  fit <- fit_initial(
    illness_death,
    models = illness.death.models,
    tau = 1
  )

  est <- tmle_estimate(fit, risk(death), scale_intensity(z = 0.5))

  expect_s3_class(est, "tmle_estimate")
  expect_named(est$estimate, c("target", "estimate", "se", "lower", "upper"))
  expect_equal(est$estimate$target, "death")
  expect_true(est$estimate$estimate > 0 && est$estimate$estimate < 1)
  expect_equal(
    c(est$estimate$lower, est$estimate$upper),
    est$estimate$estimate + c(-1, 1) * qnorm(0.975) * est$estimate$se
  )
  expect_equal(est$method, "tmle")
  expect_equal(est$target, risk(death))
  expect_equal(est$intervention, scale_intensity(z = 0.5))
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

  one_step <- tmle_estimate(
    fit,
    risk(death),
    scale_intensity(z = 0.5),
    control = list(one_step = TRUE)
  )
  expect_equal(one_step$method, "one_step")
  expect_named(one_step$estimate, names(est$estimate))
  expect_equal(one_step$initial_estimate, est$initial_estimate)
})

test_that("tmle_estimate() validates its arguments", {
  fit <- fixture.fit(30)

  expect_error(tmle_estimate(list(), risk(death)), "initial_fit")
  expect_error(tmle_estimate(fit, "death"), "risk\\(\\) or years_lost")
  expect_error(tmle_estimate(fit, risk(stroke)), "stroke.*not in initial_fit")
  expect_error(tmle_estimate(fit, risk(DEATH)), "DEATH")
  expect_error(
    tmle_estimate(fit, risk(death), list(z = 0.5)),
    "scale_intensity\\(\\) or calibrate"
  )
  expect_error(
    tmle_estimate(fit, risk(death), scale_intensity(w = 0.5)),
    "intervention is on w"
  )
  expect_error(
    tmle_estimate(fit, risk(death), a = 1),
    "a requires a treatment_model in fit_initial()",
    fixed = TRUE
  )
  expect_error(
    tmle_estimate(fit, risk(death), control = list(foo = 1)),
    "control"
  )
  expect_error(
    tmle_estimate(fit, risk(death), control = list(conv_const = 0)),
    "control\\$conv_const"
  )
  expect_error(
    tmle_estimate(fit, risk(death), control = list(truncate_weights = 0)),
    "control\\$truncate_weights"
  )
  expect_error(
    tmle_estimate(fit, risk(death), scale_intensity(z = ~ -1)),
    "non-negative"
  )
  expect_error(
    tmle_estimate(fit, risk(death), scale_intensity(z = ~undefined_var)),
    "Could not evaluate the factor of z"
  )
})

test_that("tmle_estimate() supports when, only_first and target_by_state", {
  fit <- fixture.fit(30)
  plain <- tmle_estimate(fit, risk(death), scale_intensity(z = 0.5))

  # only_first is a no-op here: death ends follow-up.
  expect_equal(
    tmle_estimate(
      fit,
      risk(death, only_first = TRUE),
      scale_intensity(z = 0.5)
    )$estimate,
    plain$estimate
  )

  before_z <- tmle_estimate(fit, risk(death, when = ~ z == 0))
  expect_lt(
    before_z$estimate$estimate,
    tmle_estimate(fit, risk(death))$estimate$estimate
  )
  expect_error(
    tmle_estimate(fit, risk(death, when = ~ L0 > 0)),
    "`when` uses L0"
  )

  expect_true(is.numeric(
    tmle_estimate(
      fit,
      risk(death),
      control = list(target_by_state = TRUE)
    )$estimate$estimate
  ))
})

test_that("tmle_estimate() truncates clever weights", {
  est <- tmle_estimate(
    fixture.fit(30),
    risk(death),
    scale_intensity(z = 3),
    control = list(truncate_weights = 1.5)
  )

  expect_lte(est$checks$weights["intervention", "q100"], 1.5 + 1e-8)
  expect_gt(est$checks$n_truncated, 0)
})

test_that("tmle_estimate() estimates under a baseline intervention 'a'", {
  fit <- fit_initial(
    fixture.data.treated(),
    models = illness.death.models,
    tau = 1,
    treatment_model = A0 ~ L0,
    arms = 0:1
  )

  est <- tmle_estimate(fit, risk(death), a = 1)
  # The treatment weights are 1 / P(A0 = 1 | L0) > 1 for the treated, and
  # separate from the censoring weights.
  expect_gt(est$checks$weights["treatment", "q100"], 1)
  expect_equal(
    est$checks$weights["censoring", ],
    tmle_estimate(fit, risk(death), a = 0)$checks$weights["censoring", ]
  )
  expect_true(est$estimate$estimate >= 0 && est$estimate$estimate <= 1)
  expect_output(print(est), "without intervention|with a = 1")

  expect_error(tmle_estimate(fit, risk(death), a = 5), "a")
})

test_that("scale_intensity() takes time-varying factors", {
  fit <- fixture.fit(30)
  constant <- tmle_estimate(fit, risk(death), scale_intensity(z = 2))

  expect_equal(
    tmle_estimate(fit, risk(death), scale_intensity(z = ~2))$estimate,
    constant$estimate
  )
  factor <- 2
  expect_equal(
    tmle_estimate(fit, risk(death), scale_intensity(z = ~factor))$estimate,
    constant$estimate
  )
  expect_equal(
    tmle_estimate(
      fit,
      risk(death),
      scale_intensity(z = ~ rep(2, length(time))),
      control = list(target_by_state = TRUE)
    )$estimate$estimate,
    tmle_estimate(
      fit,
      risk(death),
      scale_intensity(z = 2),
      control = list(target_by_state = TRUE)
    )$estimate$estimate
  )
  # Scaling by 1 after time 0.5 lies between scaling by 1 and by 2 always.
  early <- tmle_estimate(
    fit,
    risk(z),
    scale_intensity(z = ~ ifelse(time < 0.5, 2, 1))
  )
  expect_gt(
    early$estimate$estimate,
    tmle_estimate(fit, risk(z))$estimate$estimate
  )
  expect_lt(
    early$estimate$estimate,
    tmle_estimate(fit, risk(z), scale_intensity(z = 2))$estimate$estimate
  )
})

test_that("tmle_estimate() estimates the expected years lost", {
  fit <- fixture.fit(30)

  at_risk <- tmle_estimate(fit, risk(death))
  lost <- tmle_estimate(fit, years_lost(death))

  expect_equal(lost$target$estimand, "years_lost")
  # The integral of the risk up to tau = 1 is below the risk at tau.
  expect_lt(lost$estimate$estimate, at_risk$estimate$estimate)
  expect_gt(lost$estimate$estimate, 0)
  expect_output(print(lost), "estimate of years_lost\\(death\\)")
  expect_error(
    tmle_estimate(fit, years_lost(death), control = list(block_size = 0)),
    "control\\$block_size"
  )
})

test_that("tmle_estimate() only accepts a values prepared by fit_initial()", {
  fit <- fit_initial(
    fixture.data.treated(),
    models = illness.death.models,
    tau = 1,
    treatment_model = A0 ~ L0,
    arms = 1
  )

  expect_equal(tmle_estimate(fit, risk(death), a = 1)$a, 1)
  expect_error(
    tmle_estimate(fit, risk(death), a = 0),
    "a must be one of the values prepared by fit_initial(): 1.",
    fixed = TRUE
  )
  expect_error(
    tmle_estimate(
      fit,
      risk(death),
      calibrate(risk(z), relative_change(1), reference_a = 0),
      a = 1
    ),
    "reference_a must be one of"
  )
})

test_that("tmle_estimate() estimates several targets separately", {
  fit <- fixture.fit(30)
  scaled <- scale_intensity(z = 2)

  both <- tmle_estimate(fit, risk(death, z), scaled)
  death <- tmle_estimate(fit, risk(death), scaled)
  z <- tmle_estimate(fit, risk(z), scaled)

  expect_equal(both$estimate, rbind(death$estimate, z$estimate))
  expect_equal(
    both$initial_estimate,
    c(death$initial_estimate, z$initial_estimate)
  )
  expect_equal(both$eic, c(death$eic, z$eic))
  expect_named(both$checks$converged, c("death", "z"))
  expect_equal(both$checks$weights, death$checks$weights)
  expect_equal(
    coef(both),
    c(death = coef(death), z = coef(z)),
    ignore_attr = TRUE
  )
  expect_output(print(both), "z: Initial estimate")
})

test_that("scale_intensity() intervenes on several processes at once", {
  fit <- fixture.fit(30)

  # Scaling death by 1 is no intervention on it.
  expect_equal(
    tmle_estimate(
      fit,
      risk(death),
      scale_intensity(z = 0.5, death = 1)
    )$estimate,
    tmle_estimate(fit, risk(death), scale_intensity(z = 0.5))$estimate
  )
  expect_gt(
    tmle_estimate(
      fit,
      risk(death),
      scale_intensity(z = 0.5, death = 2)
    )$estimate$estimate,
    tmle_estimate(fit, risk(death), scale_intensity(z = 0.5))$estimate$estimate
  )
})
