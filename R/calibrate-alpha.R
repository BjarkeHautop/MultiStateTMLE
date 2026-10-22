#' Calibrate an intensity intervention to a target value
#'
#' Finds the intervention `alpha` on the intensity of the process `z_name`
#' such that the estimated risk of `z_name` (the exposure) equals a
#' calibration target `theta`, and estimates the risk of each `target` under
#' that calibrated intervention. The standard errors of the latter account
#' for the estimation of `alpha`.
#'
#' The calibration target `theta` is given by `value` on the `scale`:
#' * `"absolute"`: `theta = value`.
#' * `"ratio"`: `theta = value * reference`.
#' * `"difference"`: `theta = value + reference`.
#'
#' where `reference` is the exposure without intervention on `z_name`, under
#' the baseline treatment `reference_a`.
#'
#' @inheritParams tmle_alpha
#' @param z_name String. Name of the process in `initial_fit` whose
#'   intensity is intervened on and calibrated.
#' @param target Character vector. Names of the processes whose risks are
#'   estimated under the calibrated intervention.
#' @param value Number. The calibration target on the `scale`, see
#'   Description.
#' @param scale The scale of `value`: `"absolute"`, `"ratio"` or
#'   `"difference"`, see Description.
#' @param a Optional baseline treatment value (`0` or `1`) to intervene on,
#'   as in [tmle_alpha()].
#' @param reference_a Optional baseline treatment value (`0` or `1`) under
#'   which the reference exposure of a `"ratio"` or `"difference"` `scale`
#'   is estimated. Defaults to `a`; e.g. `a = 1, reference_a = 0` calibrates
#'   the exposure under treatment to that under placebo.
#' @param estimand The target parameter for `target`, as in [tmle_alpha()].
#' @param exposure_estimand The target parameter for the exposure `z_name`,
#'   that `value` and the calibration refer to, as in [tmle_alpha()].
#' @param ... Further arguments passed on to [tmle_alpha()].
#'
#' @return An object of class `calibrate_alpha`, a list with elements:
#'   * `estimate`: a `data.table` with one row per `target` and columns
#'     `target`, the risk `estimate` under the calibrated intervention, its
#'     standard error `se` (accounting for the estimation of alpha), the
#'     bounds `lower` and `upper` of a 95% Wald confidence interval, and the
#'     standard error `se_crude` treating alpha as known.
#'   * `calibration`: a `data.table` with the calibrated `alpha` and the
#'     calibration target `theta` (column `quantity`), with `estimate`, `se`,
#'     `lower` and `upper`.
#'   * `z_name`, `a`, `value`, `scale`, `reference_a`, `estimand` and
#'     `exposure_estimand`: as given.
#'   * `eic`, `eic_alpha_known`: lists, named by target, of the influence
#'     curves of the risk of `target`, accounting for the estimation of
#'     alpha and treating it as known.
#'   * `alpha_eic`: the influence curve of the calibrated alpha.
#'   * `checks`: diagnostics of the alpha search (`converged`, `dist`,
#'     `c_n`) and the estimated derivatives with respect to alpha of the risk
#'     of `z_name` (`derivative_z`) and of each target
#'     (`derivative_target`).
#'
#' @examples
#' fit <- fit_initial(
#'   illness_death,
#'   fit_types = list(
#'     z = list(
#'       delta = 2,
#'       model = ~ L0,
#'       at_risk = function(dt) dt[["z"]] == 0
#'     ),
#'     death = list(delta = 1, model = ~ L0 + z),
#'     censoring = list(delta = 0, model = ~ L0)
#'   ),
#'   tau = 1
#' )
#'
#' # Risk of death under the intervention on z that brings the risk of z
#' # down to 80% of its observed value.
#' cal <- calibrate_alpha(
#'   fit,
#'   target = "death",
#'   z_name = "z",
#'   value = 0.8,
#'   scale = "ratio"
#' )
#' cal$estimate
#' @export
calibrate_alpha <- function(
  initial_fit,
  target,
  z_name,
  value,
  scale = c("absolute", "ratio", "difference"),
  a = NULL,
  reference_a = a,
  estimand = c("risk", "years_lost"),
  exposure_estimand = c("risk", "years_lost"),
  block_size = 10,
  cores = 1,
  verbose = FALSE,
  ...
) {
  .assert_initial_fit(initial_fit)
  checkmate::assert_character(target, min.len = 1, any.missing = FALSE)
  checkmate::assert_string(z_name)
  checkmate::assert_number(value, finite = TRUE)
  scale <- match.arg(scale)
  checkmate::assert_int(a, lower = 0, upper = 1, null.ok = TRUE)
  checkmate::assert_int(reference_a, lower = 0, upper = 1, null.ok = TRUE)
  estimand <- match.arg(estimand)
  exposure_estimand <- match.arg(exposure_estimand)
  checkmate::assert_count(block_size, positive = TRUE)
  checkmate::assert_count(cores, positive = TRUE)
  checkmate::assert_flag(verbose)

  n <- initial_fit$n
  a_fixed <- a

  # Estimates of the exposure by default, as solve_alpha() needs.
  fit_at <- function(
    alpha,
    parameter = z_name,
    a = a_fixed,
    estimand = exposure_estimand
  ) {
    tmle_alpha(
      initial_fit,
      target = parameter,
      alpha = alpha,
      z_name = z_name,
      a = a,
      estimand = estimand,
      block_size = block_size,
      cores = cores,
      verbose = verbose,
      ...
    )
  }

  if (scale == "absolute") {
    theta <- value
    theta.eic <- 0
  } else {
    reference <- fit_at(alpha = 1, a = reference_a)
    if (scale == "ratio") {
      theta <- value * .point_estimate(reference)
      theta.eic <- value * .eic(reference)
    } else {
      theta <- value + .point_estimate(reference)
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

    est.deriv.auxiliary <- estimate_derivative(
      fit_at,
      alpha = est.alpha$alpha_hat,
      h = h,
      parameter = z_name
    )

    alpha.eic <- (1 / est.deriv.auxiliary) * (theta.eic - est.alpha$eic)
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
    est.deriv.auxiliary <- Inf
  }

  alpha.se <- sqrt(mean(alpha.eic^2) / n)

  out.target <- lapply(target, function(target1) {
    target.est <- fit_at(
      alpha = est.alpha$alpha_hat,
      parameter = target1,
      estimand = estimand
    )
    if (theta > 0) {
      est.deriv.target <- estimate_derivative(
        fit_at,
        alpha = est.alpha$alpha_hat,
        h = h,
        parameter = target1,
        estimand = estimand
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
  names(out.target) <- target

  per_target <- function(f) vapply(out.target, f, numeric(1))

  estimate <- .wald_table(
    target,
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

  structure(
    list(
      estimate = estimate,
      calibration = calibration,
      z_name = z_name,
      a = a,
      value = value,
      scale = scale,
      reference_a = reference_a,
      estimand = estimand,
      exposure_estimand = exposure_estimand,
      eic = lapply(out.target, function(x) x$eic),
      eic_alpha_known = lapply(out.target, function(x) .eic(x$target.est)),
      alpha_eic = alpha.eic,
      checks = list(
        converged = est.alpha$converged,
        dist = est.alpha$dist,
        c_n = est.alpha$c_n,
        derivative_z = est.deriv.auxiliary,
        derivative_target = per_target(function(x) x$est.deriv.target)
      )
    ),
    class = "calibrate_alpha"
  )
}
