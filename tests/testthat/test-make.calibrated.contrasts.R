# Cross-arm calibration decomposition from Rytgaard & van der Laan,
# supp. material Section A.1.

test_that("make.calibrated.contrasts decomposes total into mediated + calibrated", {
  set.seed(1)
  n <- 200

  treatment.eic <- rnorm(n, sd = 0.1)
  placebo.eic <- rnorm(n, sd = 0.1)
  calibrated.eic <- rnorm(n, sd = 0.1)

  treatment.fit <- list(
    estimate = c(tmle.est = 0.7),
    eic = treatment.eic
  )
  placebo.fit <- list(
    estimate = c(tmle.est = 0.3),
    eic = placebo.eic
  )
  calibrated.fit <- list(
    estimate = c(target.est = 0.5),
    eic = calibrated.eic
  )

  out <- make.calibrated.contrasts(
    treatment.fit = treatment.fit,
    placebo.fit = placebo.fit,
    calibrated.fit = calibrated.fit
  )

  expect_equal(out$contrast, c("total", "mediated", "calibrated"))
  expect_equal(out$estimate[out$contrast == "total"], 0.7 - 0.3)
  expect_equal(out$estimate[out$contrast == "mediated"], 0.7 - 0.5)
  expect_equal(out$estimate[out$contrast == "calibrated"], 0.5 - 0.3)

  expect_equal(
    out$estimate[out$contrast == "total"],
    out$estimate[out$contrast == "mediated"] +
      out$estimate[out$contrast == "calibrated"]
  )
  expect_equal(attr(out, "decomposition.error"), 0)
  expect_equal(attr(out, "eic.decomposition.error"), 0)
})

test_that("make.calibrated.contrasts computes se/CI from the influence curves", {
  n <- 50
  treatment.eic <- rep(c(1, -1), length.out = n)
  placebo.eic <- rep(0, n)
  calibrated.eic <- rep(c(0.5, -0.5, 0.2, -0.2), length.out = n)

  treatment.fit <- list(estimate = c(tmle.est = 1), eic = treatment.eic)
  placebo.fit <- list(estimate = c(tmle.est = 0), eic = placebo.eic)
  calibrated.fit <- list(estimate = c(target.est = 0.4), eic = calibrated.eic)

  out <- make.calibrated.contrasts(
    treatment.fit = treatment.fit,
    placebo.fit = placebo.fit,
    calibrated.fit = calibrated.fit,
    conf.level = 0.95
  )

  crit <- qnorm(0.975)

  expected.se <- c(
    total = sd(treatment.eic - placebo.eic),
    mediated = sd(treatment.eic - calibrated.eic),
    calibrated = sd(calibrated.eic - placebo.eic)
  ) / sqrt(n)

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

test_that("make.calibrated.contrasts attributes the total correctly at the decomposition's endpoints", {
  n <- 30
  treatment.eic <- rnorm(n, sd = 0.1)
  placebo.eic <- rnorm(n, sd = 0.1)

  treatment.fit <- list(estimate = c(tmle.est = 0.8), eic = treatment.eic)
  placebo.fit <- list(estimate = c(tmle.est = 0.2), eic = placebo.eic)

  calibrated.fit.at.treatment <- list(
    estimate = c(target.est = 0.8),
    eic = treatment.eic
  )
  out.at.treatment <- make.calibrated.contrasts(
    treatment.fit = treatment.fit,
    placebo.fit = placebo.fit,
    calibrated.fit = calibrated.fit.at.treatment
  )
  expect_equal(out.at.treatment$estimate[out.at.treatment$contrast == "mediated"], 0)
  expect_equal(
    out.at.treatment$estimate[out.at.treatment$contrast == "calibrated"],
    out.at.treatment$estimate[out.at.treatment$contrast == "total"]
  )

  calibrated.fit.at.placebo <- list(
    estimate = c(target.est = 0.2),
    eic = placebo.eic
  )
  out.at.placebo <- make.calibrated.contrasts(
    treatment.fit = treatment.fit,
    placebo.fit = placebo.fit,
    calibrated.fit = calibrated.fit.at.placebo
  )
  expect_equal(out.at.placebo$estimate[out.at.placebo$contrast == "calibrated"], 0)
  expect_equal(
    out.at.placebo$estimate[out.at.placebo$contrast == "mediated"],
    out.at.placebo$estimate[out.at.placebo$contrast == "total"]
  )
})

test_that("make.calibrated.contrasts errors on mismatched eic lengths", {
  treatment.fit <- list(estimate = c(tmle.est = 1), eic = rep(0, 10))
  placebo.fit <- list(estimate = c(tmle.est = 0), eic = rep(0, 5))
  calibrated.fit <- list(estimate = c(target.est = 0), eic = rep(0, 10))

  expect_error(
    make.calibrated.contrasts(
      treatment.fit = treatment.fit,
      placebo.fit = placebo.fit,
      calibrated.fit = calibrated.fit
    ),
    "must have the same length"
  )
})

test_that("make.calibrated.contrasts supports selecting a named target", {
  treatment.fit <- list(estimate = c(tmle.est = 1), eic = rep(0.1, 20))
  placebo.fit <- list(estimate = c(tmle.est = 0), eic = rep(-0.1, 20))
  calibrated.fit <- list(
    estimate = c("target.est.arm1" = 0.4, "target.est.arm2" = 0.6),
    eic = list(arm1 = rep(0, 20), arm2 = rep(0.05, 20))
  )

  out <- make.calibrated.contrasts(
    target = "arm1",
    treatment.fit = treatment.fit,
    placebo.fit = placebo.fit,
    calibrated.fit = calibrated.fit
  )

  expect_equal(out$estimate[out$contrast == "mediated"], 1 - 0.4)
  expect_equal(out$estimate[out$contrast == "calibrated"], 0.4 - 0)
  expect_true(all(out$target == "arm1"))
})
