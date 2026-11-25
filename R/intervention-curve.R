#' Trace the intervention curve over a grid of intensity scalings
#'
#' For each scaling `alpha` of the intensity of the exposure process,
#' estimates the exposure and the target under [scale_intensity()] with
#' that `alpha`. Plotting the target estimates against the exposure
#' estimates traces out how the target responds to changes in the exposure
#' on its own scale, such as its risk. The standard errors of the target
#' estimates account for calibrating to the estimated exposure rather than
#' to `alpha`, as in [calibrate()].
#'
#' @inheritParams tmle_estimate
#' @param target The target parameter: [risk()] or [years_lost()] of one or
#'   more processes.
#' @param exposure The exposure: [risk()] or [years_lost()] of the single
#'   process whose intensity is scaled.
#' @param alpha Non-negative numeric vector. The scalings to evaluate.
#'
#' @return An object of class `intervention_curve`, a list with elements:
#'   * `estimate`: a `data.table` with one row per value of `alpha` and
#'     target process, and columns `alpha`, `target`, the exposure
#'     (`exposure_estimate`, with standard error `exposure_se`), and the
#'     target (`estimate`, with standard error `se` and the bounds `lower`
#'     and `upper` of a 95% Wald confidence interval).
#'   * `target`, `exposure` and `a`: as given.
#'   * `checks`: a `data.table` with the estimated derivatives with respect
#'     to `alpha` of the exposure (`derivative_exposure`) and target
#'     (`derivative_target`) estimates, per `alpha` and target process.
#'   * `eic`: a list with, per value of `alpha`, the influence curves of the
#'     exposure estimate (`exposure_eic`), of `alpha` calibrated to it
#'     (`alpha_eic`), and lists named by target process of the target
#'     estimates accounting for its calibration (`eic`) and treating `alpha`
#'     as known (`eic_alpha_known`), as in [tmle_estimate()].
#'
#' @examples
#' fit <- fit_initial(
#'   illness_death,
#'   models = list(
#'     hazard_model(illness ~ biomarker + treatment, limit = 1),
#'     death ~ biomarker + treatment + illness,
#'     censoring ~ biomarker
#'   ),
#'   tau = 1
#' )
#'
#' \donttest{
#' # How the risk of death responds to the risk of illness, as the intensity
#' # of illness is scaled from 0 (no illness) to 2 (twice the observed
#' # intensity):
#' curve <- intervention_curve(
#'   fit,
#'   risk(death),
#'   exposure = risk(illness),
#'   alpha = c(0, 0.5, 1, 2)
#' )
#' curve
#'
#' # Plot the risk of death against the risk of illness, with 95% confidence
#' # intervals:
#' est <- curve$estimate
#' plot(
#'   est$exposure_estimate,
#'   est$estimate,
#'   type = "b",
#'   ylim = range(est$lower, est$upper),
#'   xlab = "Risk of illness by time 1",
#'   ylab = "Risk of death by time 1"
#' )
#' arrows(
#'   est$exposure_estimate,
#'   est$lower,
#'   est$exposure_estimate,
#'   est$upper,
#'   angle = 90,
#'   code = 3,
#'   length = 0.03
#' )
#' }
#' @export
intervention_curve <- function(
  initial_fit,
  target,
  exposure,
  alpha = seq(0, 5, length.out = 10),
  a = NULL,
  control = list(),
  verbose = FALSE
) {
  .assert_initial_fit(initial_fit)
  .assert_target(target, initial_fit)
  .assert_target(exposure, initial_fit, "exposure")
  if (length(exposure$processes) != 1) {
    stop("`exposure` must be of a single process.")
  }
  checkmate::assert_numeric(
    alpha,
    lower = 0,
    finite = TRUE,
    any.missing = FALSE,
    min.len = 1
  )
  .assert_a(a, initial_fit)
  control <- .tmle_control(control)
  checkmate::assert_flag(verbose)

  n <- initial_fit$n
  z_name <- exposure$processes
  alpha_grid <- alpha

  fit_at <- function(alpha, parameter) {
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

  alpha.out <- lapply(alpha_grid, function(alpha) {
    if (alpha > 0) {
      exposure.est <- fit_at(alpha = alpha, parameter = exposure)

      h <- 0.3 * n^(-1 / 6) * sqrt(mean(.eic(exposure.est)^2))

      est.deriv.exposure <- estimate_derivative(
        fit_at,
        alpha = alpha,
        h = h,
        parameter = exposure
      )

      alpha.eic <- (1 / est.deriv.exposure) * (0 - .eic(exposure.est))
    } else {
      exposure.est <- list(
        estimate = data.table(estimate = 0, se = 0),
        eic = list(0)
      )
      alpha.eic <- 0
      est.deriv.exposure <- Inf
    }

    out.target <- lapply(target$processes, function(process) {
      parameter <- target
      parameter$processes <- process
      target.est <- fit_at(alpha = alpha, parameter = parameter)

      if (alpha > 0) {
        est.deriv.target <- estimate_derivative(
          fit_at,
          alpha = alpha,
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

    per_target <- function(f) vapply(out.target, f, numeric(1))

    estimate <- cbind(
      data.table(
        alpha = alpha,
        exposure_estimate = .point_estimate(exposure.est),
        exposure_se = exposure.est$estimate$se[1]
      ),
      .wald_table(
        target$processes,
        per_target(function(x) .point_estimate(x$target.est)),
        per_target(function(x) sqrt(mean(x$eic^2) / n))
      )
    )
    setcolorder(estimate, c("alpha", "target"))

    checks <- data.table(
      alpha = alpha,
      target = target$processes,
      derivative_exposure = est.deriv.exposure,
      derivative_target = per_target(function(x) x$est.deriv.target)
    )

    eic <- list(
      exposure_eic = .eic(exposure.est),
      alpha_eic = alpha.eic,
      eic = lapply(out.target, function(x) x$eic),
      eic_alpha_known = lapply(out.target, function(x) .eic(x$target.est))
    )

    list(estimate = estimate, checks = checks, eic = eic)
  })

  structure(
    list(
      estimate = rbindlist(lapply(alpha.out, function(x) x$estimate)),
      target = target,
      exposure = exposure,
      a = a,
      checks = rbindlist(lapply(alpha.out, function(x) x$checks)),
      eic = stats::setNames(
        lapply(alpha.out, function(x) x$eic),
        paste0("alpha = ", alpha_grid)
      )
    ),
    class = "intervention_curve"
  )
}
