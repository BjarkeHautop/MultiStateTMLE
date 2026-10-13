#' Calibrate an intensity intervention to a target value
#'
#' Finds the intervention `alpha` on the intensity of the process `z_name`
#' such that the estimated risk of `z_name` equals a calibration target, and
#' estimates the risk of `target` under that calibrated intervention. The
#' standard error of the latter accounts for the estimation of `alpha`.
#'
#' Exactly one of `theta`, `rho`, `delta` and `rho_1a` specifies the
#' calibration target:
#' * `theta`: the risk of `z_name` itself.
#' * `rho`: a fraction of the risk of `z_name` without intervention.
#' * `delta`: an absolute difference from the risk of `z_name` without
#'   intervention.
#' * `rho_1a`: a fraction of the risk of `z_name` without intervention on
#'   `z_name`, under the opposite baseline treatment `1 - a`.
#'
#' @inheritParams tmle_alpha
#' @param target Character vector. Names of the processes whose risks are
#'   estimated under the calibrated intervention.
#' @param theta,rho,delta,rho_1a Numbers specifying the calibration target,
#'   see Description. Give exactly one.
#' @param output_eic Logical. If `TRUE`, include the influence curves in the
#'   output.
#' @param ... Further arguments passed on to [tmle_alpha()].
#'
#' @return A list with elements:
#'   * `estimate`: a named vector with the calibrated `alpha.est` and its
#'     `alpha.se`, the risk of `target` (`target.est`) with standard error
#'     `target.se` (accounting for the estimation of alpha) and
#'     `target.se.crude` (treating alpha as known), and the calibration
#'     target `theta` with standard error `theta.se`.
#'   * `checks`: diagnostics of the alpha search and derivative estimates.
#'   * `alpha.eic`, `target.eic` and `eic` (if `output_eic = TRUE`): the
#'     influence curves of the calibrated alpha, of the risk of `target`
#'     treating alpha as known, and of the risk of `target`.
#'
#'   With several targets, the `target.*` entries of `estimate` and `checks`
#'   are suffixed by the target name (e.g. `target.est.outcome1`), and
#'   `target.eic` and `eic` are lists named by target.
#'
#' @examples
#' fit <- fit_initial(
#'   illness_death,
#'   fit_types = list(
#'     z = list(
#'       model = "Surv(tstart, tstop, delta == 2) ~ L0",
#'       at_risk = function(dt) dt[["z"]] == 0
#'     ),
#'     outcome1 = "Surv(tstart, tstop, delta == 1) ~ L0 + z",
#'     censoring = "Surv(tstart, tstop, delta == 0) ~ L0"
#'   ),
#'   tau = 1
#' )
#'
#' # Risk of outcome1 under the intervention on z that brings the risk of z
#' # down to 80% of its observed value.
#' cal <- calibrate_alpha(fit, target = "outcome1", z_name = "z", rho = 0.8)
#' cal$estimate
#' @export
calibrate_alpha <- function(
  initial_fit,
  target,
  z_name = "z",
  theta = NULL,
  rho = NULL,
  delta = NULL,
  rho_1a = NULL,
  a = NULL,
  output_eic = FALSE,
  cores = 1,
  verbose = FALSE,
  ...
) {
  .assert_initial_fit(initial_fit)
  checkmate::assert_character(target, min.len = 1, any.missing = FALSE)
  checkmate::assert_string(z_name)
  checkmate::assert_number(theta, lower = 0, finite = TRUE, null.ok = TRUE)
  checkmate::assert_number(rho, lower = 0, finite = TRUE, null.ok = TRUE)
  checkmate::assert_number(delta, finite = TRUE, null.ok = TRUE)
  checkmate::assert_number(rho_1a, lower = 0, finite = TRUE, null.ok = TRUE)
  n_targets <- length(c(theta, rho, delta, rho_1a))
  if (n_targets != 1) {
    stop("Give exactly one of theta, rho, delta and rho_1a.")
  }
  checkmate::assert_int(a, lower = 0, upper = 1, null.ok = TRUE)
  if (!is.null(rho_1a) && is.null(a)) {
    stop("rho_1a requires a baseline treatment a.")
  }
  checkmate::assert_flag(output_eic)
  checkmate::assert_count(cores, positive = TRUE)
  checkmate::assert_flag(verbose)

  n <- length(unique(initial_fit[["tmp.long"]][["id"]]))
  a_fixed <- a

  fit_at <- function(alpha, parameter = z_name, a = a_fixed) {
    tmle_alpha(
      initial_fit,
      target = parameter,
      alpha = alpha,
      z_name = z_name,
      a = a,
      output_eic = TRUE,
      cores = cores,
      verbose = verbose,
      ...
    )
  }

  if (!is.null(rho)) {
    est.aux.1 <- fit_at(alpha = 1)
    theta <- rho * .point_estimate(est.aux.1)
    theta.eic <- rho * est.aux.1$eic
  } else if (!is.null(rho_1a)) {
    est.aux.1 <- fit_at(alpha = 1, a = 1 - a)
    theta <- rho_1a * .point_estimate(est.aux.1)
    theta.eic <- rho_1a * est.aux.1$eic
  } else if (!is.null(delta)) {
    est.aux.1 <- fit_at(alpha = 1)
    theta <- delta + .point_estimate(est.aux.1)
    theta.eic <- est.aux.1$eic
  } else {
    theta.eic <- 0
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

    est.deriv.auxiliary <- estimate_derivative(
      fit_at,
      alpha = est.alpha$alpha.hat,
      h = h,
      parameter = z_name
    )

    alpha.eic <- (1 / est.deriv.auxiliary) * (theta.eic - est.alpha$eic)
  } else {
    est.alpha <- list(
      alpha.hat = 0,
      eic = 0,
      converged = TRUE,
      dist = 0,
      c.n = 0,
      grid = data.frame(alpha = 0, psi = 0)
    )
    alpha.eic <- 0
    est.deriv.auxiliary <- Inf
  }

  alpha.se <- sqrt(mean(alpha.eic^2) / n)

  out.target <- lapply(target, function(target1) {
    target.est <- fit_at(alpha = est.alpha$alpha.hat, parameter = target1)
    if (theta > 0) {
      est.deriv.target <- estimate_derivative(
        fit_at,
        alpha = est.alpha$alpha.hat,
        h = h,
        parameter = target1
      )
      eic <- target.est[["eic"]] + est.deriv.target * alpha.eic
    } else {
      est.deriv.target <- Inf
      eic <- target.est[["eic"]]
    }
    list(
      target.est = target.est,
      est.deriv.target = est.deriv.target,
      eic = eic
    )
  })
  names(out.target) <- if (length(target) > 1) target else ""

  per_target <- function(f) vapply(out.target, f, numeric(1))

  out <- list(
    estimate = c(
      alpha.est = est.alpha$alpha.hat,
      alpha.se = alpha.se,
      target.est = per_target(function(x) .point_estimate(x$target.est)),
      target.se = per_target(function(x) sqrt(mean(x$eic^2) / n)),
      target.se.crude = per_target(function(x) {
        as.numeric(x$target.est[["estimate"]][["se"]])
      }),
      theta = theta,
      theta.se = theta.se
    ),
    checks = c(
      est.deriv.auxiliary = est.deriv.auxiliary,
      est.alpha.converged = est.alpha$converged,
      est.alpha.dist = est.alpha$dist,
      est.alpha.cn = est.alpha$c.n,
      est.deriv.target = per_target(function(x) x$est.deriv.target)
    )
  )

  if (output_eic) {
    target.eic <- lapply(out.target, function(x) x$target.est[["eic"]])
    eic <- lapply(out.target, function(x) x$eic)
    if (length(target) == 1) {
      target.eic <- target.eic[[1]]
      eic <- eic[[1]]
    }
    out$alpha.eic <- alpha.eic
    out$target.eic <- target.eic
    out$eic <- eic
  }

  out
}
