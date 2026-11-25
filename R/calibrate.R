# tmle_estimate() under a calibrate() intervention: solve for the scaling
# alpha of the exposure process at which the exposure attains theta, given
# by the calibration rule, and estimate the targets under it, with
# influence curves accounting for alpha being estimated.
.tmle_calibrated <- function(
  initial_fit,
  target,
  intervention,
  a,
  control,
  verbose
) {
  exposure <- intervention$exposure
  z_name <- exposure$processes
  rule <- intervention$toward
  reference_a <- if (is.null(intervention$reference_a)) {
    a
  } else {
    intervention$reference_a
  }
  n <- initial_fit$n
  a_fixed <- a

  # The fit of `parameter` (a tmle_target, the exposure by default, as
  # solve_alpha() needs) for its first process, under scaling alpha of illness.
  fit_at <- function(alpha, parameter = exposure, a = a_fixed) {
    .tmle_fixed(
      initial_fit,
      parameter$processes[1],
      parameter,
      .scale_one(z_name, alpha),
      a,
      control,
      verbose
    )
  }

  if (rule$rule == "level") {
    theta <- rule$value
    theta.eic <- 0
  } else {
    reference <- fit_at(alpha = 1, a = reference_a)
    if (rule$rule == "relative_change") {
      theta <- rule$value * .point_estimate(reference)
      theta.eic <- rule$value * .eic(reference)
    } else {
      theta <- rule$value + .point_estimate(reference)
      theta.eic <- .eic(reference)
    }
  }

  theta.se <- sqrt(mean(theta.eic^2) / n)

  if (theta < 0) {
    stop("Inadmissible calibration target: theta = ", theta, ".")
  }

  if (theta > 0) {
    est.alpha <- solve_alpha(
      fun = fit_at,
      theta = theta,
      c_n = n^(-1 / 2) / log(n),
      verbose = verbose
    )

    if (!isTRUE(est.alpha$converged)) {
      stop(
        "Calibration failed to converge: no alpha achieves theta = ",
        theta,
        " (closest achieved: ",
        est.alpha$grid$psi[which.min(abs(est.alpha$grid$psi - theta))],
        ", dist = ",
        est.alpha$dist,
        "); theta may exceed L(P), the achievable ceiling of Psi_z^alpha(P)."
      )
    }

    h <- 0.3 * n^(-1 / 6) * sqrt(mean(est.alpha$eic^2))

    est.deriv.exposure <- estimate_derivative(
      fit_at,
      alpha = est.alpha$alpha_hat,
      h = h
    )

    alpha.eic <- (1 / est.deriv.exposure) * (theta.eic - est.alpha$eic)
  } else {
    est.alpha <- list(
      alpha_hat = 0,
      eic = 0,
      converged = TRUE,
      dist = 0,
      c_n = 0,
      grid = data.frame(alpha = 0, psi = 0)
    )
    alpha.eic <- 0
    est.deriv.exposure <- Inf
  }

  alpha.se <- sqrt(mean(alpha.eic^2) / n)

  out.target <- lapply(target$processes, function(process) {
    parameter <- target
    parameter$processes <- process
    target.est <- fit_at(alpha = est.alpha$alpha_hat, parameter = parameter)
    if (theta > 0) {
      est.deriv.target <- estimate_derivative(
        fit_at,
        alpha = est.alpha$alpha_hat,
        h = h,
        parameter = parameter
      )
      eic <- .eic(target.est) + est.deriv.target * alpha.eic
    } else {
      est.deriv.target <- Inf
      eic <- .eic(target.est)
    }
    list(
      target.est = target.est,
      est.deriv.target = est.deriv.target,
      eic = eic
    )
  })
  names(out.target) <- target$processes
  fits <- lapply(out.target, function(x) x$target.est)

  per_target <- function(f) vapply(out.target, f, numeric(1))

  estimate <- .wald_table(
    target$processes,
    per_target(function(x) .point_estimate(x$target.est)),
    per_target(function(x) sqrt(mean(x$eic^2) / n))
  )
  set(
    estimate,
    j = "se_crude",
    value = per_target(function(x) x$target.est[["estimate"]][["se"]])
  )

  calibration <- .wald_table(
    c("alpha", "theta"),
    c(est.alpha$alpha_hat, theta),
    c(alpha.se, theta.se)
  )
  setnames(calibration, "target", "quantity")

  checks <- .combine_targets(unname(fits))$checks
  checks$calibration <- list(
    converged = est.alpha$converged,
    dist = est.alpha$dist,
    c_n = est.alpha$c_n,
    derivative_exposure = est.deriv.exposure,
    derivative_target = per_target(function(x) x$est.deriv.target)
  )

  structure(
    list(
      estimate = estimate,
      calibration = calibration,
      target = target,
      intervention = intervention,
      a = a,
      method = fits[[1]]$method,
      eic = lapply(out.target, function(x) x$eic),
      eic_alpha_known = lapply(out.target, function(x) .eic(x$target.est)),
      alpha_eic = alpha.eic,
      checks = checks
    ),
    class = "tmle_estimate"
  )
}
