# Recurrent processes: a process that can jump any number of times, so that
# risk() is the expected number of events, the state space is capped at
# max_count, and risk(only_first = TRUE) differs from the default.

test_that("fit_initial() caps the state space of a recurrent process at max_count", {
  skip_if_not_installed("simevent")

  fit <- recurrent.fit(80, max_count = 3)
  expect_equal(fit$processes$type, c("recurrent", "terminal", "terminal"))
  expect_equal(fit$depend_matrix$relapse, 0:3)
  # The last state jumps to itself.
  expect_equal(fit$depend_matrix$s.relapse, c(2L, 3L, 4L, 4L))
  # Each row is in the state of its number of relapses so far, capped.
  tl <- fit$tmp_long
  expect_equal(tl$state, pmin(tl$relapse, 3) + 1)

  capped <- recurrent.fit(80, max_count = 1)
  expect_equal(capped$depend_matrix$relapse, 0:1)
  expect_equal(capped$depend_matrix$s.relapse, c(2L, 2L))
  tl <- capped$tmp_long
  expect_true(any(tl$relapse > 1))
  expect_equal(tl$state, pmin(tl$relapse, 1) + 1)
})

test_that(".compute_q_clever() counts recurrent events beyond the capped state", {
  # States Z = 0, 1, where state 2 (Z >= 1) jumps to itself.
  states <- data.table::data.table(
    state = c(1L, 2L),
    Z = c(0, 1),
    s.Z = c(2L, 2L)
  )
  dt_id <- data.table::data.table(
    state = c(1L, 1L),
    state.row.index = c(1L, 1L),
    P.Z.1 = c(0.2, 0.2),
    P.Z.2 = c(0.3, 0.3),
    P.D.1 = c(0.1, 0.1),
    P.D.2 = c(0.1, 0.1)
  )

  out <- .compute_q_clever(
    dt_id = dt_id,
    states = states,
    process_types = list(Z = "recurrent", D = "terminal"),
    parameter = "Z"
  )

  # Q_2(s) is the expected number of Z jumps at the last time: 0.2, 0.3.
  # Q_1(1) = 0.2 + 0.7 * 0.2 + 0.2 * 0.3 = 0.40
  # Q_1(2) = 0.3 + 0.6 * 0.3 + 0.3 * 0.3 = 0.57, the jump from the capped
  #   state counting and staying in it.
  expect_equal(out$Q.1, c(0.40, 0.2))
  expect_equal(out$Q.2, c(0.57, 0.3))
  expect_equal(out$Q, c(0.40, 0.2))
  # A jump of the recurrent target adds 1 to the successor's Q.
  expect_equal(out$clever.Q.Z1, c(0.3 + 1, 0.3 + 1))
  expect_equal(out$clever.Q.Z0, c(0.2, 0.2))
})

test_that("risk() of a recurrent process is the expected number of events", {
  skip_if_not_installed("simevent")
  fit <- recurrent.fit(80)

  count <- tmle_estimate(fit, risk(relapse))
  first <- tmle_estimate(fit, risk(relapse, only_first = TRUE))
  expect_true(all(c(count$checks$converged, first$checks$converged)))
  expect_gt(first$estimate$estimate, 0)
  expect_lt(first$estimate$estimate, 1)
  # E[N(tau)] > P(N(tau) >= 1) once some individuals relapse twice.
  expect_gt(count$estimate$estimate, first$estimate$estimate)

  # With max_count = 1 the state space no longer tells relapses apart, but
  # every relapse still counts.
  capped <- recurrent.fit(80, max_count = 1)
  expect_gt(
    tmle_estimate(capped, risk(relapse))$estimate$estimate,
    tmle_estimate(capped, risk(relapse, only_first = TRUE))$estimate$estimate
  )

  # Several targets are each estimated as on their own.
  both <- tmle_estimate(
    fit,
    risk(death, relapse),
    scale_intensity(relapse = 0.5)
  )
  expect_equal(
    both$estimate,
    rbind(
      tmle_estimate(fit, risk(death), scale_intensity(relapse = 0.5))$estimate,
      tmle_estimate(fit, risk(relapse), scale_intensity(relapse = 0.5))$estimate
    )
  )
})

test_that("risk(only_first = TRUE) is the expected number of relapses from the state without any", {
  skip_if_not_installed("simevent")
  # Two routes to the risk of a first relapse: as a one-jump target, and as
  # the expected number of relapses counted only when none came before.
  # They agree in the plug-in, the targeted estimate and the efficient
  # influence curve, and so need the Q recursion of every individual to
  # run up to tau, also after their second relapse.
  fit <- recurrent.fit(80)
  for (intervention in list(NULL, scale_intensity(relapse = 0.5))) {
    first <- tmle_estimate(
      fit,
      risk(relapse, only_first = TRUE),
      intervention
    )
    from_none <- tmle_estimate(
      fit,
      risk(relapse, when = ~ relapse == 0),
      intervention
    )
    expect_equal(first$initial_estimate, from_none$initial_estimate)
    expect_equal(first$estimate, from_none$estimate)
    expect_equal(first$eic, from_none$eic)
  }
})

test_that("scale_intensity() of a recurrent process moves its targets monotonically", {
  skip_if_not_installed("simevent")
  fit <- recurrent.fit(80)

  # The estimates of target under the factor factor(alpha) on relapse, for
  # increasing alpha.
  estimates <- function(target, factor = identity, initial = FALSE) {
    vapply(
      c(0.5, 1, 2),
      function(alpha) {
        est <- tmle_estimate(
          fit,
          target,
          scale_intensity(relapse = factor(alpha))
        )
        if (initial) est$initial_estimate[[1]] else est$estimate$estimate
      },
      numeric(1)
    )
  }
  expect_true(all(diff(estimates(risk(relapse))) > 0))
  expect_true(all(diff(estimates(risk(relapse, only_first = TRUE))) > 0))
  first_only <- function(alpha) ~ ifelse(relapse == 0, alpha, 1)
  expect_true(all(diff(estimates(risk(relapse), first_only)) > 0))
  # Relapses raise the fitted death hazard, so the plug-in risk of death
  # grows with alpha. (The targeted estimates need not, with n = 80.)
  expect_true(all(diff(estimates(risk(death), initial = TRUE)) > 0))

  # Removing relapses altogether leaves none.
  none <- scale_intensity(relapse = 0)
  expect_equal(tmle_estimate(fit, risk(relapse), none)$estimate$estimate, 0)
  expect_equal(
    tmle_estimate(
      fit,
      risk(relapse, only_first = TRUE),
      none
    )$estimate$estimate,
    0
  )
  # Then every death is a death before relapse.
  expect_equal(
    tmle_estimate(fit, risk(death, when = ~ relapse == 0), none)$estimate,
    tmle_estimate(fit, risk(death), none)$estimate
  )
})

test_that("only_first and interventions need a recurrent process in the state space", {
  skip_if_not_installed("simevent")
  # No hazard depends on relapse, so its count is not tracked.
  fit <- fit_initial(
    simulate.recurrent(80, seed = 1),
    models = list(relapse ~ L0, death ~ L0, censoring ~ L0),
    tau = 1
  )
  expect_equal(fit$processes$type[1], "recurrent")

  # Its expected number of events needs no state space.
  expect_gt(tmle_estimate(fit, risk(relapse))$estimate$estimate, 0)
  expect_error(
    tmle_estimate(fit, risk(relapse, only_first = TRUE)),
    "only_first needs the number of events of relapse"
  )
  expect_error(
    tmle_estimate(fit, risk(death), scale_intensity(relapse = 0.5)),
    "To intervene on relapse"
  )
})

test_that("a factor that depends on the count leaves later relapses as observed", {
  skip_if_not_installed("simevent")
  fit <- recurrent.fit(80)

  initial <- function(factor) {
    tmle_estimate(
      fit,
      risk(relapse),
      scale_intensity(relapse = factor)
    )$initial_estimate
  }
  # Halving the first relapse only: later relapses are at their observed
  # rate, above the halved one.
  first_only <- initial(~ ifelse(relapse == 0, 0.5, 1))
  expect_gt(first_only, initial(0.5))
  expect_lt(first_only, initial(1))
})

test_that("a factor that depends on the count applies in each state", {
  skip_if_not_installed("simevent")
  fit <- recurrent.fit(80)

  # Only relapses from state relapse == 0 matter for the first relapse, and
  # there the factor is 0.5 either way.
  constant <- tmle_estimate(
    fit,
    risk(relapse, only_first = TRUE),
    scale_intensity(relapse = 0.5)
  )
  by_count <- tmle_estimate(
    fit,
    risk(relapse, only_first = TRUE),
    scale_intensity(relapse = ~ ifelse(relapse == 0, 0.5, 1))
  )
  expect_equal(by_count$initial_estimate, constant$initial_estimate)
  expect_equal(by_count$estimate, constant$estimate)
})

test_that("a factor sees the counts capped at max_count", {
  skip_if_not_installed("simevent")
  # With max_count = 1, relapse >= 1 and relapse == 1 are the same state,
  # also on the observed rows of individuals with more relapses.
  capped <- recurrent.fit(80, max_count = 1)
  at_least <- tmle_estimate(
    capped,
    risk(relapse),
    scale_intensity(relapse = ~ ifelse(relapse >= 1, 0.5, 1))
  )
  exactly <- tmle_estimate(
    capped,
    risk(relapse),
    scale_intensity(relapse = ~ ifelse(relapse == 1, 0.5, 1))
  )
  expect_equal(exactly$estimate, at_least$estimate)
})

test_that(".update_alpha_weights() accumulates the likelihood ratio of a time-varying factor", {
  # One individual who relapses at time 0.3 under a factor of 0.5, which
  # is 1 from time 0.6 on.
  tl <- data.table::data.table(
    id = 1,
    time = c(0.1, 0.3, 0.6, 0.8),
    delta = c(0L, 1L, 0L, 0L),
    relapse = c(0L, 0L, 1L, 1L),
    at.risk.relapse = 1,
    P.relapse = 0.1,
    alpha.relapse = c(0.5, 0.5, 1, 1)
  )
  .update_alpha_weights(tl, "relapse")

  # The weight at each time is the product over the earlier times of
  # alpha^dN * exp(-(alpha - 1) * dLambda): it keeps the factor 0.5 of the
  # relapse and the cumulative hazard up to then after alpha returns to 1.
  ratio <- 0.5^(tl$delta == 1) * exp(-(0.5 - 1) * 0.1)
  expect_equal(
    tl$clever.weight.alpha,
    c(1, ratio[1], ratio[1] * ratio[2], ratio[1] * ratio[2])
  )
})

# Simulation-based recovery of the closed-form truth of helper-recurrent.R.
test_that("tmle_estimate() recovers recurrent targets", {
  skip_on_cran()
  skip_if_not_installed("simevent")

  fit <- recurrent.fit(300, seed = 1)
  truth <- list(
    `1` = recurrent.truth.constant(1),
    `0.5` = recurrent.truth.constant(0.5),
    `2` = recurrent.truth.constant(2)
  )
  cases <- list(
    list(risk(relapse), 1, "relapse"),
    list(risk(relapse), 0.5, "relapse"),
    list(risk(relapse), 2, "relapse"),
    list(risk(relapse, only_first = TRUE), 1, "relapse_first"),
    list(risk(relapse, only_first = TRUE), 0.5, "relapse_first"),
    list(risk(death), 0.5, "death"),
    list(risk(death), 2, "death"),
    list(risk(death, when = ~ relapse == 0), 0.5, "death_before_relapse"),
    list(years_lost(death), 0.5, "years_lost_death")
  )
  for (case in cases) {
    target <- case[[1]]
    alpha <- case[[2]]
    est <- tmle_estimate(
      fit,
      target,
      scale_intensity(relapse = alpha)
    )$estimate
    expect_lt(
      abs(est$estimate - truth[[as.character(alpha)]][[case[[3]]]]) / est$se,
      4,
      label = paste(format(target), "at alpha =", alpha)
    )
  }
})
