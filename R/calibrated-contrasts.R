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
#'   as returned by [tmle_alpha()] with different baseline treatments `a`.
#' @param calibrated_fit Fit for the calibrated intervention, as returned by
#'   [calibrate_alpha()] with the `a` of `treatment_fit`. It may contain
#'   estimates for several targets; the one of `treatment_fit` is used.
#' @param level Number between 0 and 1. Confidence level of the intervals.
#'
#' @return An object of class `calibrated_contrasts`, a list with elements:
#'   * `estimate`: a `data.table` with one row per contrast (`total`,
#'     `mediated` and `calibrated`) and columns `target`, `contrast`,
#'     `estimate`, `se`, `lower` and `upper`.
#'   * `target` and `level`: the target and the confidence level.
#'   * `checks`: a list with `decomposition_error` and
#'     `eic_decomposition_error`, how far the estimates and influence curves
#'     are from satisfying total = mediated + calibrated, which should be
#'     zero up to numerical error.
#'
#' @examples
#' # illness_death with a randomised treatment A0 added, for illustration.
#' set.seed(1)
#' treated <- illness_death
#' treated$A0 <- rbinom(max(treated$id), 1, 0.5)[treated$id]
#'
#' fit <- fit_initial(
#'   treated,
#'   fit_types = list(
#'     z = list(
#'       delta = 2,
#'       model = ~ L0 + A0,
#'       at_risk = function(dt) dt[["z"]] == 0
#'     ),
#'     death = list(delta = 1, model = ~ L0 + A0 + z),
#'     censoring = list(delta = 0, model = ~L0)
#'   ),
#'   tau = 1,
#'   treatment_model = A0 ~ L0
#' )
#'
#' # How much of the effect of treatment on death is mediated by its effect
#' # on the risk of z?
#' calibrated_contrasts(
#'   treatment_fit = tmle_alpha(fit, target = "death", a = 1),
#'   placebo_fit = tmle_alpha(fit, target = "death", a = 0),
#'   calibrated_fit = calibrate_alpha(
#'     fit,
#'     target = "death",
#'     z_name = "z",
#'     value = 1,
#'     scale = "ratio",
#'     a = 1,
#'     reference_a = 0
#'   )
#' )
#' @export
calibrated_contrasts <- function(
  treatment_fit,
  placebo_fit,
  calibrated_fit,
  level = 0.95
) {
  checkmate::assert_class(treatment_fit, "tmle_alpha")
  checkmate::assert_class(placebo_fit, "tmle_alpha")
  checkmate::assert_class(calibrated_fit, "calibrate_alpha")
  checkmate::assert_number(level, lower = 0, upper = 1)

  target <- treatment_fit$estimate$target
  if (!identical(placebo_fit$estimate$target, target)) {
    stop("treatment_fit and placebo_fit must have the same target.")
  }
  if (
    is.null(treatment_fit$a) ||
      is.null(placebo_fit$a) ||
      treatment_fit$a == placebo_fit$a
  ) {
    stop(
      "treatment_fit and placebo_fit must be under different baseline ",
      "treatments a."
    )
  }
  if (!identical(calibrated_fit$a, treatment_fit$a)) {
    stop(
      "calibrated_fit must be under the baseline treatment of ",
      "treatment_fit, a = ",
      treatment_fit$a,
      "."
    )
  }
  estimands <- c(
    treatment_fit$estimand,
    placebo_fit$estimand,
    calibrated_fit$estimand
  )
  if (length(unique(estimands)) != 1) {
    stop(
      "treatment_fit, placebo_fit and calibrated_fit must have the same ",
      "estimand."
    )
  }
  calibrated.row <- match(target, calibrated_fit$estimate$target)
  if (is.na(calibrated.row)) {
    stop(
      "The target of treatment_fit, '",
      target,
      "', is not among the targets of calibrated_fit: ",
      paste(calibrated_fit$estimate$target, collapse = ", "),
      "."
    )
  }

  treatment.est <- .point_estimate(treatment_fit)
  placebo.est <- .point_estimate(placebo_fit)
  calibrated.est <- calibrated_fit$estimate$estimate[calibrated.row]

  treatment.eic <- .eic(treatment_fit)
  placebo.eic <- .eic(placebo_fit)
  calibrated.eic <- calibrated_fit$eic[[calibrated.row]]

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

    se <- sqrt(mean(eic^2 / n))

    critical.value <-
      stats::qnorm(1 - (1 - level) / 2)

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

  structure(
    list(
      estimate = out,
      target = target,
      level = level,
      checks = list(
        decomposition_error = decomposition.error,
        eic_decomposition_error = eic.decomposition.error
      )
    ),
    class = "calibrated_contrasts"
  )
}
