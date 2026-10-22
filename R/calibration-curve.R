#' Trace the calibration curve over a grid of interventions
#'
#' For each intervention `alpha` in `alpha_grid` on the intensity of the
#' process `z_name`, estimates the risk of `z_name` (the exposure) and of
#' each `target` under that intervention. Plotting the target estimates
#' against the exposure estimates traces out the calibration curve relating
#' the two. The standard errors of the target estimates account for
#' calibrating to the estimated exposure risk rather than to `alpha`.
#'
#' @inheritParams tmle_alpha
#' @param z_name String. Name of the process in `initial_fit` whose
#'   intensity `alpha_grid` scales.
#' @param target Character vector. Names of the processes whose risks are
#'   estimated.
#' @param alpha_grid Non-negative numeric vector. The interventions to
#'   evaluate.
#' @param min_iter Non-negative integer. Minimum number of TMLE iterations,
#'   passed on to [tmle_alpha()].
#' @param block_size_z,block_size_target Optional positive integers. If
#'   given, the expected time lost is estimated for `z_name` and `target`,
#'   respectively, instead of the risk; see `years_lost` in [tmle_alpha()].
#' @param ... Further arguments passed on to [tmle_alpha()].
#'
#' @return A list with elements:
#'   * `estimate`: a `data.table` with one row per value of `alpha_grid` and
#'     columns `alpha`, `exposure_est`, `exposure_se`, `target_est` and
#'     `target_se` (suffixed by `_<target>` when there are several
#'     targets).
#'   * `checks`: per value of `alpha_grid`, the estimated derivatives of the
#'     exposure and target estimates with respect to `alpha`.
#'   * `eics`: per value of `alpha_grid`, the influence curves `alpha_eic`,
#'     `calibration_eic`, `target_eic` and `eic`.
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
#' curve <- calibration_curve(
#'   fit,
#'   target = "death",
#'   z_name = "z",
#'   alpha_grid = c(0, 0.5, 1)
#' )
#' curve$estimate
#' @export
calibration_curve <- function(
  initial_fit,
  target,
  z_name,
  alpha_grid = seq(0, 5, length.out = 10),
  a = NULL,
  min_iter = 1,
  target_by_state = FALSE,
  block_size_z = NULL,
  block_size_target = NULL,
  cores = 1,
  verbose = FALSE,
  ...
) {
  .assert_initial_fit(initial_fit)
  checkmate::assert_character(target, min.len = 1, any.missing = FALSE)
  checkmate::assert_string(z_name)
  checkmate::assert_numeric(
    alpha_grid,
    lower = 0,
    finite = TRUE,
    any.missing = FALSE,
    min.len = 1
  )
  checkmate::assert_int(a, lower = 0, upper = 1, null.ok = TRUE)
  checkmate::assert_count(min_iter)
  checkmate::assert_flag(target_by_state)
  checkmate::assert_count(block_size_z, positive = TRUE, null.ok = TRUE)
  checkmate::assert_count(block_size_target, positive = TRUE, null.ok = TRUE)
  checkmate::assert_count(cores, positive = TRUE)
  checkmate::assert_flag(verbose)

  n <- length(unique(initial_fit[["tmp_long"]][["id"]]))

  fit_at <- function(alpha, parameter = z_name, block_size = NULL) {
    tmle_alpha(
      initial_fit,
      target = parameter,
      alpha = alpha,
      z_name = z_name,
      a = a,
      years_lost = block_size,
      min_iter = min_iter,
      target_by_state = target_by_state,
      cores = cores,
      verbose = verbose,
      ...
    )
  }

  alpha.out <- lapply(alpha_grid, function(alpha) {
    if (alpha > 0) {
      calibration.est <- fit_at(
        alpha = alpha,
        parameter = z_name,
        block_size = block_size_z
      )

      h <- 0.3 * n^(-1 / 6) * sqrt(mean(calibration.est$eic^2))

      est.deriv.auxiliary <- estimate_derivative(
        fit_at,
        alpha = alpha,
        h = h,
        parameter = z_name,
        block_size = block_size_z
      )

      alpha.eic <- (1 / est.deriv.auxiliary) * (0 - calibration.est$eic)
    } else {
      calibration.est <- list(estimate = c(tmle_est = 0, se = 0), eic = 0)
      alpha.eic <- 0
      est.deriv.auxiliary <- Inf
    }

    out.target <- lapply(target, function(target1) {
      target.est <- fit_at(
        alpha = alpha,
        parameter = target1,
        block_size = block_size_target
      )

      if (alpha > 0) {
        est.deriv.target <- estimate_derivative(
          fit_at,
          alpha = alpha,
          h = h,
          parameter = target1,
          block_size = block_size_target
        )

        target.se <- sqrt(
          mean((target.est[["eic"]] + est.deriv.target * alpha.eic)^2) / n
        )
      } else {
        est.deriv.target <- Inf
        target.se <- sqrt(mean((target.est[["eic"]])^2) / n)
      }

      list(
        target.est = target.est,
        target.se = target.se,
        est.deriv.target = est.deriv.target
      )
    })

    names(out.target) <- if (length(target) > 1) target else ""

    out.estimate <- c(
      alpha = alpha,
      exposure_est = .point_estimate(calibration.est),
      exposure_se = calibration.est$estimate[["se"]],
      .target_names(
        vapply(
          out.target,
          function(x) .point_estimate(x[["target.est"]]),
          numeric(1)
        ),
        "target_est"
      ),
      .target_names(
        vapply(out.target, function(x) x[["target.se"]], numeric(1)),
        "target_se"
      )
    )

    out.checks <- c(
      est_deriv_auxiliary = est.deriv.auxiliary,
      .target_names(
        vapply(out.target, function(x) x[["est.deriv.target"]], numeric(1)),
        "est_deriv_target"
      )
    )

    target.eic <- lapply(out.target, function(x) x[["target.est"]][["eic"]])
    eic <- lapply(out.target, function(x) {
      if (alpha > 0) {
        x[["target.est"]][["eic"]] + x[["est.deriv.target"]] * alpha.eic
      } else {
        x[["target.est"]][["eic"]]
      }
    })
    if (length(target) == 1) {
      target.eic <- target.eic[[1]]
      eic <- eic[[1]]
    }
    list(
      estimate = out.estimate,
      checks = out.checks,
      alpha_eic = alpha.eic,
      calibration_eic = calibration.est$eic,
      target_eic = target.eic,
      eic = eic
    )
  })

  names(alpha.out) <- paste0("alpha = ", alpha_grid)

  estimate.table <- rbindlist(
    lapply(alpha.out, function(x) as.list(x$estimate)),
    fill = TRUE
  )

  checks.list <- lapply(alpha.out, function(x) x$checks)

  eic.list <- lapply(alpha.out, function(x) {
    x[c("alpha_eic", "calibration_eic", "target_eic", "eic")]
  })

  list(estimate = estimate.table, checks = checks.list, eics = eic.list)
}
