# Cross-arm calibration decomposition from Rytgaard & van der Laan,
# supp. material Section A.1.

test_that("calibrated_contrasts() decomposes total into mediated + calibrated", {
  set.seed(1)
  n <- 200

  treatment.eic <- rnorm(n, sd = 0.1)
  placebo.eic <- rnorm(n, sd = 0.1)
  calibrated.eic <- rnorm(n, sd = 0.1)

  treatment.fit <- fake.tmle.fit(0.7, treatment.eic)
  placebo.fit <- fake.tmle.fit(0.3, placebo.eic, a = 0)
  calibrated.fit <- fake.calibrated.fit(0.5, calibrated.eic)

  res <- calibrated_contrasts(
    treatment_fit = treatment.fit,
    placebo_fit = placebo.fit,
    calibrated_fit = calibrated.fit
  )
  expect_s3_class(res, "calibrated_contrasts")
  expect_output(print(res), "Decomposition of the effect of a on death")
  out <- res$estimate

  expect_named(
    out,
    c("target", "contrast", "estimate", "se", "lower", "upper")
  )
  expect_true(all(out$target == "death"))
  expect_equal(out$contrast, c("total", "mediated", "calibrated"))
  expect_equal(out$estimate[out$contrast == "total"], 0.7 - 0.3)
  expect_equal(out$estimate[out$contrast == "mediated"], 0.7 - 0.5)
  expect_equal(out$estimate[out$contrast == "calibrated"], 0.5 - 0.3)

  expect_equal(
    out$estimate[out$contrast == "total"],
    out$estimate[out$contrast == "mediated"] +
      out$estimate[out$contrast == "calibrated"]
  )
  expect_equal(res$checks$decomposition_error, 0)
  expect_equal(res$checks$eic_decomposition_error, 0)
})

test_that("calibrated_contrasts() computes se/CI from the influence curves", {
  n <- 50
  treatment.eic <- rep(c(1, -1), length.out = n)
  placebo.eic <- rep(0, n)
  calibrated.eic <- rep(c(0.5, -0.5, 0.2, -0.2), length.out = n)

  treatment.fit <- fake.tmle.fit(1, treatment.eic)
  placebo.fit <- fake.tmle.fit(0, placebo.eic, a = 0)
  calibrated.fit <- fake.calibrated.fit(0.4, calibrated.eic)

  out <- calibrated_contrasts(
    treatment_fit = treatment.fit,
    placebo_fit = placebo.fit,
    calibrated_fit = calibrated.fit,
    level = 0.95
  )$estimate

  crit <- qnorm(0.975)

  expected.se <- c(
    total = sqrt(mean((treatment.eic - placebo.eic)^2 / n)),
    mediated = sqrt(mean((treatment.eic - calibrated.eic)^2 / n)),
    calibrated = sqrt(mean((calibrated.eic - placebo.eic)^2 / n))
  )

  for (row in c("total", "mediated", "calibrated")) {
    this.row <- out[out$contrast == row]
    expect_equal(this.row$se, unname(expected.se[row]))
    expect_equal(
      this.row$lower,
      this.row$estimate - crit * expected.se[row],
      ignore_attr = TRUE
    )
    expect_equal(
      this.row$upper,
      this.row$estimate + crit * expected.se[row],
      ignore_attr = TRUE
    )
  }
})

test_that("calibrated_contrasts() attributes the total correctly at the decomposition's endpoints", {
  n <- 30
  treatment.eic <- rnorm(n, sd = 0.1)
  placebo.eic <- rnorm(n, sd = 0.1)

  treatment.fit <- fake.tmle.fit(0.8, treatment.eic)
  placebo.fit <- fake.tmle.fit(0.2, placebo.eic, a = 0)

  calibrated.fit.at.treatment <- fake.calibrated.fit(0.8, treatment.eic)
  out.at.treatment <- calibrated_contrasts(
    treatment_fit = treatment.fit,
    placebo_fit = placebo.fit,
    calibrated_fit = calibrated.fit.at.treatment
  )$estimate
  expect_equal(
    out.at.treatment$estimate[out.at.treatment$contrast == "mediated"],
    0
  )
  expect_equal(
    out.at.treatment$estimate[out.at.treatment$contrast == "calibrated"],
    out.at.treatment$estimate[out.at.treatment$contrast == "total"]
  )

  calibrated.fit.at.placebo <- fake.calibrated.fit(0.2, placebo.eic)
  out.at.placebo <- calibrated_contrasts(
    treatment_fit = treatment.fit,
    placebo_fit = placebo.fit,
    calibrated_fit = calibrated.fit.at.placebo
  )$estimate
  expect_equal(
    out.at.placebo$estimate[out.at.placebo$contrast == "calibrated"],
    0
  )
  expect_equal(
    out.at.placebo$estimate[out.at.placebo$contrast == "mediated"],
    out.at.placebo$estimate[out.at.placebo$contrast == "total"]
  )
})

test_that("calibrated_contrasts() errors on mismatched eic lengths", {
  treatment.fit <- fake.tmle.fit(1, rep(0, 10))
  placebo.fit <- fake.tmle.fit(0, rep(0, 5), a = 0)
  calibrated.fit <- fake.calibrated.fit(0, rep(0, 10))

  expect_error(
    calibrated_contrasts(
      treatment_fit = treatment.fit,
      placebo_fit = placebo.fit,
      calibrated_fit = calibrated.fit
    ),
    "must have the same length"
  )
})

test_that("calibrated_contrasts() checks that the fits belong together", {
  eic <- c(0.1, -0.1)
  treatment.fit <- fake.tmle.fit(0.7, eic)
  placebo.fit <- fake.tmle.fit(0.3, eic, a = 0)
  calibrated.fit <- fake.calibrated.fit(0.5, eic)

  expect_error(
    calibrated_contrasts(
      treatment.fit,
      fake.tmle.fit(0.3, eic, target = "z", a = 0),
      calibrated.fit
    ),
    "same target"
  )
  expect_error(
    calibrated_contrasts(
      treatment.fit,
      fake.tmle.fit(0.3, eic),
      calibrated.fit
    ),
    "different baseline treatments"
  )
  expect_error(
    calibrated_contrasts(
      treatment.fit,
      fake.tmle.fit(0.3, eic, a = NULL),
      calibrated.fit
    ),
    "different baseline treatments"
  )
  expect_error(
    calibrated_contrasts(
      treatment.fit,
      placebo.fit,
      fake.calibrated.fit(0.5, eic, a = 0)
    ),
    "baseline treatment of treatment_fit, a = 1"
  )
  expect_error(
    calibrated_contrasts(
      treatment.fit,
      placebo.fit,
      fake.calibrated.fit(0.5, eic, estimand = "years_lost")
    ),
    "same estimand"
  )
})

test_that("calibrated_contrasts() defaults to the treatment fit's target", {
  eic <- c(0.1, -0.1)
  calibrated.fit <- fake.calibrated.fit(
    c(0.5, 0.9),
    list(eic, eic),
    target = c("z", "death")
  )

  out <- calibrated_contrasts(
    fake.tmle.fit(0.7, eic),
    fake.tmle.fit(0.3, eic, a = 0),
    calibrated.fit
  )$estimate
  expect_equal(out$target, rep("death", 3))
  expect_equal(out$estimate, c(0.7 - 0.3, 0.7 - 0.9, 0.9 - 0.3))

  expect_error(
    calibrated_contrasts(
      fake.tmle.fit(0.7, eic, target = "stroke"),
      fake.tmle.fit(0.3, eic, target = "stroke", a = 0),
      calibrated.fit
    ),
    "not among the targets"
  )
  expect_error(
    calibrated_contrasts(list(), fake.tmle.fit(0.3, eic), calibrated.fit),
    "tmle_alpha"
  )
})
