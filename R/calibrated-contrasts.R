#' Decompose a treatment effect into mediated and calibrated contrasts
#'
#' Given estimates of the risk of an outcome under treatment, under placebo,
#' and under a calibrated intervention (as returned by [tmle_alpha()] and
#' [calibrate_alpha()]), computes the total treatment effect and its
#' decomposition into a mediated contrast (treatment versus calibrated) and
#' a calibrated contrast (calibrated versus placebo), with Wald confidence
#' intervals based on the influence curves.
#'
#' @param treatment_fit,placebo_fit Fits for the treatment and placebo arms,
#'   as returned by [tmle_alpha()] with `output_eic = TRUE`.
#' @param calibrated_fit Fit for the calibrated intervention, as returned by
#'   [calibrate_alpha()] with `output_eic = TRUE`.
#' @param target Optional string. If `calibrated_fit` contains estimates for
#'   several targets, the name of the one to use.
#' @param conf_level Number between 0 and 1. Confidence level of the
#'   intervals.
#'
#' @return A `data.table` with one row per contrast (`total`, `mediated`
#'   and `calibrated`) and columns `estimate`, `se`, `lower` and `upper`.
#'   The attributes `decomposition.error` and `eic.decomposition.error` give
#'   how far the estimates and influence curves are from satisfying
#'   total = mediated + calibrated, which should be zero up to numerical
#'   error.
#'
#' @examples
#' set.seed(1405)
#' treatment_fit <- list(estimate = c(tmle.est = 0.7), eic = rnorm(200, sd = 0.1))
#' placebo_fit <- list(estimate = c(tmle.est = 0.3), eic = rnorm(200, sd = 0.1))
#' calibrated_fit <- list(estimate = c(target.est = 0.5), eic = rnorm(200, sd = 0.1))
#'
#' calibrated_contrasts(treatment_fit, placebo_fit, calibrated_fit)
#' @export
calibrated_contrasts <- function(
  treatment_fit,
  placebo_fit,
  calibrated_fit,
  target = NULL,
  conf_level = 0.95
) {
  checkmate::assert_list(treatment_fit)
  checkmate::assert_list(placebo_fit)
  checkmate::assert_list(calibrated_fit)
  checkmate::assert_string(target, null.ok = TRUE)
  checkmate::assert_number(conf_level, lower = 0, upper = 1)

  treatment.est <-
    treatment_fit$estimate[["tmle.est"]]

  placebo.est <-
    placebo_fit$estimate[["tmle.est"]]

  calibrated.est <-
    calibrated_fit$estimate[[ifelse(
      length(target) > 0,
      paste0("target.est.", target),
      "target.est"
    )]]

  treatment.eic <- treatment_fit$eic
  placebo.eic <- placebo_fit$eic
  calibrated.eic <- if (length(target) > 0) {
    calibrated_fit$eic[[target]]
  } else {
    calibrated_fit$eic
  }

  n <- length(treatment.eic)

  if (
    length(placebo.eic) != n ||
      length(calibrated.eic) != n
  ) {
    stop("The EIC vectors must have the same length.")
  }

  contrast <- function(
    estimate.1,
    estimate.0,
    eic.1,
    eic.0
  ) {
    estimate <- estimate.1 - estimate.0
    eic <- eic.1 - eic.0

    se <- stats::sd(eic) / sqrt(n)

    critical.value <-
      stats::qnorm(1 - (1 - conf_level) / 2)

    c(
      estimate = estimate,
      se = se,
      lower = estimate - critical.value * se,
      upper = estimate + critical.value * se
    )
  }

  total <- contrast(
    estimate.1 = treatment.est,
    estimate.0 = placebo.est,
    eic.1 = treatment.eic,
    eic.0 = placebo.eic
  )

  mediated <- contrast(
    estimate.1 = treatment.est,
    estimate.0 = calibrated.est,
    eic.1 = treatment.eic,
    eic.0 = calibrated.eic
  )

  calibrated <- contrast(
    estimate.1 = calibrated.est,
    estimate.0 = placebo.est,
    eic.1 = calibrated.eic,
    eic.0 = placebo.eic
  )

  out <- rbind(total, mediated, calibrated)
  out <- data.table(target = target, contrast = rownames(out), out)

  decomposition.error <-
    total[["estimate"]] -
    mediated[["estimate"]] -
    calibrated[["estimate"]]

  eic.decomposition.error <-
    max(abs(
      (treatment.eic - placebo.eic) -
        (treatment.eic - calibrated.eic) -
        (calibrated.eic - placebo.eic)
    ))

  attr(out, "decomposition.error") <-
    decomposition.error

  attr(out, "eic.decomposition.error") <-
    eic.decomposition.error

  out
}
