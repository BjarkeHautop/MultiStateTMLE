#' Stochastic intensity interventions
#'
#' Describe the intervention that [tmle_estimate()] estimates the target
#' under.
#'
#' `scale_intensity()` is a fixed intensity-scaling intervention: it
#' multiplies the intensity of each named process by its factor `alpha`,
#' \eqn{\Lambda^z \mapsto \alpha \Lambda^z}. `alpha = 1` is no intervention,
#' `alpha < 1` attenuates the process and `alpha > 1` intensifies it.
#'
#' `calibrate()` is a calibrated intervention: the intensity of the process
#' of `exposure` is scaled by the `alpha` that brings the exposure to the
#' value given by `toward`. The target is then a composite parameter, whose
#' standard error accounts for `alpha` being estimated. The calibration
#' rules are:
#' * `level(theta)`: the exposure equals `theta`.
#' * `absolute_change(delta)`: the exposure is `delta` plus its value
#'   without intervention.
#' * `relative_change(rho)`: the exposure is `rho` times its value without
#'   intervention.
#'
#' @param ... Named factors, one per process to intervene on, such as
#'   `z = 0.5`. Each is a non-negative number, or a one-sided formula giving
#'   a time-varying factor, such as `z = ~ ifelse(time < 0.5, 0.5, 1)`, that
#'   may use `time`, baseline covariates and the processes (their number of
#'   events so far, where counts above the `max_count` of [fit_initial()]
#'   count as `max_count`).
#' @param exposure The exposure to calibrate, a [risk()] or [years_lost()]
#'   of the single process whose intensity is scaled.
#' @param toward The calibration rule: `level()`, `absolute_change()` or
#'   `relative_change()`.
#' @param reference_a Optional baseline treatment value (`0` or `1`) under
#'   which the exposure without intervention is estimated, for
#'   `absolute_change()` and `relative_change()`. Defaults to the `a` given
#'   to [tmle_estimate()]. A different value calibrates across arms, e.g.
#'   the risk of the exposure under treatment to that under placebo.
#' @param theta Number between 0 and 1. The exposure to attain.
#' @param delta Number. The change in the exposure.
#' @param rho Non-negative number. The factor on the exposure.
#'
#' @return `scale_intensity()` and `calibrate()` return an object of class
#'   `tmle_intervention`; `level()`, `absolute_change()` and
#'   `relative_change()` one of class `calibration_rule`.
#'
#' @examples
#' # Fixed intensity scalings. Halve the intensity of z, remove z altogether,
#' # or double it:
#' scale_intensity(z = 0.5)
#' scale_intensity(z = 0)
#' scale_intensity(z = 2)
#'
#' # A time-varying scaling: halve the intensity of z before time 0.5 only.
#' scale_intensity(z = ~ ifelse(time < 0.5, 0.5, 1))
#'
#' # A scaling that depends on a baseline covariate: halve it only for
#' # individuals with L0 > 0.
#' scale_intensity(z = ~ ifelse(L0 > 0, 0.5, 1))
#'
#' # A scaling that depends on the number of events so far: halve the
#' # intensity of a recurrent process until its first event only.
#' scale_intensity(relapse = ~ ifelse(relapse == 0, 0.5, 1))
#'
#' # Several processes at once:
#' scale_intensity(z = 0.5, death = 0.8)
#'
#' # Calibrated interventions, scaling the intensity of z so that its risk
#' # equals 0.1, falls by 0.05, or falls to 80% of its value without
#' # intervention:
#' calibrate(risk(z), level(0.1))
#' calibrate(risk(z), absolute_change(-0.05))
#' calibrate(risk(z), relative_change(0.8))
#'
#' # Use them in tmle_estimate():
#' fit <- fit_initial(
#'   illness_death,
#'   models = list(
#'     hazard_model(z ~ L0, limit = 1),
#'     death ~ L0 + z,
#'     censoring ~ L0
#'   ),
#'   tau = 1
#' )
#' tmle_estimate(fit, risk(death), scale_intensity(z = ~ ifelse(L0 > 0, 0.5, 1)))
#' \donttest{
#' tmle_estimate(fit, risk(death), calibrate(risk(z), level(0.1)))
#' }
#' @name interventions
NULL

#' @rdname interventions
#' @export
scale_intensity <- function(...) {
  alpha <- list(...)
  if (length(alpha) == 0) {
    stop("scale_intensity() needs a factor for at least one process.")
  }
  if (is.null(names(alpha)) || any(names(alpha) == "")) {
    stop(
      "The factors in scale_intensity() must be named by process, ",
      "e.g. scale_intensity(z = 0.5)."
    )
  }
  if (anyDuplicated(names(alpha))) {
    stop("Each process can only be given once to scale_intensity().")
  }
  for (name in names(alpha)) {
    factor <- alpha[[name]]
    if (inherits(factor, "formula")) {
      if (length(factor) != 2) {
        stop("The factor of ", name, " must be a one-sided formula.")
      }
    } else {
      checkmate::assert_number(
        factor,
        lower = 0,
        finite = TRUE,
        .var.name = name
      )
    }
  }
  structure(
    list(type = "scale_intensity", alpha = alpha),
    class = "tmle_intervention"
  )
}

#' @rdname interventions
#' @export
calibrate <- function(exposure, toward, reference_a = NULL) {
  if (!inherits(exposure, "tmle_target")) {
    stop("`exposure` must be given by risk() or years_lost().")
  }
  if (length(exposure$processes) != 1) {
    stop("`exposure` must be of a single process.")
  }
  if (!inherits(toward, "calibration_rule")) {
    stop(
      "`toward` must be given by level(), absolute_change() or ",
      "relative_change()."
    )
  }
  checkmate::assert_choice(reference_a, 0:1, null.ok = TRUE)
  if (toward$rule == "level" && !is.null(reference_a)) {
    stop(
      "`reference_a` is only used by absolute_change() and relative_change()."
    )
  }
  structure(
    list(
      type = "calibrate",
      exposure = exposure,
      toward = toward,
      reference_a = reference_a
    ),
    class = "tmle_intervention"
  )
}

#' @rdname interventions
#' @export
level <- function(theta) {
  checkmate::assert_number(theta, lower = 0, upper = 1)
  .calibration_rule("level", theta)
}

#' @rdname interventions
#' @export
absolute_change <- function(delta) {
  checkmate::assert_number(delta, finite = TRUE)
  .calibration_rule("absolute_change", delta)
}

#' @rdname interventions
#' @export
relative_change <- function(rho) {
  checkmate::assert_number(rho, lower = 0, finite = TRUE)
  .calibration_rule("relative_change", rho)
}

.calibration_rule <- function(rule, value) {
  structure(list(rule = rule, value = value), class = "calibration_rule")
}

#' @export
format.calibration_rule <- function(x, ...) {
  paste0(x$rule, "(", format(x$value), ")")
}

#' @export
format.tmle_intervention <- function(x, ...) {
  if (x$type == "scale_intensity") {
    factors <- vapply(
      x$alpha,
      function(factor) {
        if (inherits(factor, "formula")) {
          .deparse_one_line(factor)
        } else {
          format(factor)
        }
      },
      character(1)
    )
    args <- paste(names(x$alpha), "=", factors)
    return(paste0("scale_intensity(", paste(args, collapse = ", "), ")"))
  }
  args <- c(format(x$exposure), format(x$toward))
  if (!is.null(x$reference_a)) {
    args <- c(args, paste("reference_a =", x$reference_a))
  }
  paste0("calibrate(", paste(args, collapse = ", "), ")")
}

#' @export
print.tmle_intervention <- function(x, ...) {
  cat("Intervention: ", format(x), "\n", sep = "")
  invisible(x)
}

#' @export
print.calibration_rule <- function(x, ...) {
  cat("Calibration rule: ", format(x), "\n", sep = "")
  invisible(x)
}

# Check that the processes an intervention refers to are in `initial_fit`.
.assert_intervention <- function(intervention, initial_fit) {
  if (is.null(intervention)) {
    return(invisible(NULL))
  }
  if (!inherits(intervention, "tmle_intervention")) {
    stop(
      "`intervention` must be given by scale_intensity() or calibrate()."
    )
  }
  processes <- if (intervention$type == "scale_intensity") {
    names(intervention$alpha)
  } else {
    intervention$exposure$processes
  }
  unknown <- setdiff(processes, initial_fit$process_names)
  if (length(unknown) > 0) {
    stop(
      "The intervention is on ",
      paste(unknown, collapse = ", "),
      ", which is not a process in initial_fit: ",
      paste(initial_fit$process_names, collapse = ", "),
      "."
    )
  }
  if (intervention$type == "calibrate") {
    .assert_a(intervention$reference_a, initial_fit, "reference_a")
  }
  invisible(intervention)
}
