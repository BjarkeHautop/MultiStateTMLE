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
#' @param estimand The target parameter for `target`, as in [tmle_alpha()].
#' @param exposure_estimand The target parameter for the exposure `z_name`,
#'   as in [tmle_alpha()].
#' @param control Settings for the estimates, passed on to [tmle_alpha()].
#'
#' @return An object of class `calibration_curve`, a list with elements:
#'   * `estimate`: a `data.table` with one row per value of `alpha_grid` and
#'     `target`, and columns `alpha`, `target`, the risk of `z_name`
#'     (`exposure_estimate`, with standard error `exposure_se`), and the risk
#'     of `target` (`estimate`, with standard error `se` and the bounds
#'     `lower` and `upper` of a 95% Wald confidence interval).
#'   * `z_name`, `a`, `estimand` and `exposure_estimand`: as given.
#'   * `checks`: a `data.table` with the estimated derivatives with respect
#'     to `alpha` of the exposure (`derivative_exposure`) and target
#'     (`derivative_target`) estimates, per `alpha` and `target`.
#'   * `eic`: a list with, per value of `alpha_grid`, the influence curves of
#'     the exposure estimate (`exposure_eic`), of `alpha` calibrated to it
#'     (`alpha_eic`), and lists named by target of the target estimates
#'     accounting for its calibration (`eic`) and treating `alpha` as known
#'     (`eic_alpha_known`), as in [calibrate_alpha()].
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
  estimand = c("risk", "years_lost"),
  exposure_estimand = c("risk", "years_lost"),
  control = list(),
  cores = 1,
  verbose = FALSE
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
  .assert_a(a, initial_fit)
  estimand <- match.arg(estimand)
  exposure_estimand <- match.arg(exposure_estimand)
  .tmle_control(control)
  checkmate::assert_count(cores, positive = TRUE)
  checkmate::assert_flag(verbose)

  n <- initial_fit$n

  fit_at <- function(alpha, parameter, estimand) {
    tmle_alpha(
      initial_fit,
      target = parameter,
      alpha = alpha,
      z_name = z_name,
      a = a,
      estimand = estimand,
      control = control,
      cores = cores,
      verbose = verbose
    )
  }

  alpha.out <- lapply(alpha_grid, function(alpha) {
    if (alpha > 0) {
      calibration.est <- fit_at(
        alpha = alpha,
        parameter = z_name,
        estimand = exposure_estimand
      )

      h <- 0.3 * n^(-1 / 6) * sqrt(mean(.eic(calibration.est)^2))

      est.deriv.auxiliary <- estimate_derivative(
        fit_at,
        alpha = alpha,
        h = h,
        parameter = z_name,
        estimand = exposure_estimand
      )

      alpha.eic <- (1 / est.deriv.auxiliary) * (0 - .eic(calibration.est))
    } else {
      calibration.est <- list(
        estimate = data.table(estimate = 0, se = 0),
        eic = list(0)
      )
      alpha.eic <- 0
      est.deriv.auxiliary <- Inf
    }

    out.target <- lapply(target, function(target1) {
      target.est <- fit_at(
        alpha = alpha,
        parameter = target1,
        estimand = estimand
      )

      if (alpha > 0) {
        est.deriv.target <- estimate_derivative(
          fit_at,
          alpha = alpha,
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

    estimate <- cbind(
      data.table(
        alpha = alpha,
        exposure_estimate = .point_estimate(calibration.est),
        exposure_se = calibration.est$estimate$se[1]
      ),
      .wald_table(
        target,
        per_target(function(x) .point_estimate(x$target.est)),
        per_target(function(x) sqrt(mean(x$eic^2) / n))
      )
    )
    setcolorder(estimate, c("alpha", "target"))

    checks <- data.table(
      alpha = alpha,
      target = target,
      derivative_exposure = est.deriv.auxiliary,
      derivative_target = per_target(function(x) x$est.deriv.target)
    )

    eic <- list(
      exposure_eic = .eic(calibration.est),
      alpha_eic = alpha.eic,
      eic = lapply(out.target, function(x) x$eic),
      eic_alpha_known = lapply(out.target, function(x) .eic(x$target.est))
    )

    list(estimate = estimate, checks = checks, eic = eic)
  })

  structure(
    list(
      estimate = rbindlist(lapply(alpha.out, function(x) x$estimate)),
      z_name = z_name,
      a = a,
      estimand = estimand,
      exposure_estimand = exposure_estimand,
      checks = rbindlist(lapply(alpha.out, function(x) x$checks)),
      eic = stats::setNames(
        lapply(alpha.out, function(x) x$eic),
        paste0("alpha = ", alpha_grid)
      )
    ),
    class = "calibration_curve"
  )
}
