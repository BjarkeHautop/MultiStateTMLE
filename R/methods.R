#' Methods for the result objects
#'
#' Print, extract the estimates of, and compute confidence intervals for the
#' objects returned by [fit_initial()], [tmle_alpha()], [calibrate_alpha()],
#' [calibration_curve()] and [calibrated_contrasts()].
#'
#' @param x,object A result object.
#' @param parm Optional character vector of targets to compute confidence
#'   intervals for. Defaults to all.
#' @param level Number between 0 and 1. Confidence level.
#' @param ... Ignored.
#'
#' @return `print()` returns `x` invisibly. `coef()` returns a numeric vector
#'   of estimates named by target, and `confint()` a matrix of Wald
#'   confidence intervals with one row per target.
#'
#' @examples
#' fit <- fit_initial(
#'   illness_death,
#'   fit_types = list(
#'     z = list(
#'       delta = 2,
#'       model = ~L0,
#'       at_risk = function(dt) dt[["z"]] == 0
#'     ),
#'     death = list(delta = 1, model = ~ L0 + z),
#'     censoring = list(delta = 0, model = ~L0)
#'   ),
#'   tau = 1
#' )
#' fit
#'
#' est <- tmle_alpha(fit, target = "death", alpha = 0.5, z_name = "z")
#' est
#' coef(est)
#' confint(est, level = 0.9)
#' @name result-methods
NULL

#' @rdname result-methods
#' @export
print.initial_fit <- function(x, ...) {
  cat(
    "Initial fit of ",
    x$n,
    " individuals up to tau = ",
    x$tau,
    "\n\n",
    sep = ""
  )
  print(x$processes)
  if (!is.null(x$treatment_model)) {
    cat(
      "\nTreatment model: ",
      .deparse_one_line(x$treatment_model),
      ", prepared for arms a = ",
      paste(x$arms, collapse = ", "),
      "\n",
      sep = ""
    )
  }
  invisible(x)
}

#' @rdname result-methods
#' @export
print.tmle_alpha <- function(x, ...) {
  method <- if (x$method == "one_step") "One-step" else "TMLE"
  cat(
    method,
    " estimate ",
    if (identical(x$estimand, "years_lost")) "of years lost ",
    .describe_intervention(x),
    "\n\n",
    sep = ""
  )
  print(x$estimate)
  cat("\n")
  targets <- names(x$initial_estimate)
  for (target in targets) {
    if (length(targets) > 1) {
      cat(target, ": ", sep = "")
    }
    cat("Initial estimate:", format(x$initial_estimate[[target]], digits = 4))
    if (x$method == "tmle") {
      iter <- x$checks$iter[[target]]
      if (x$checks$converged[[target]]) {
        cat("; converged after", iter, "iterations")
      } else {
        cat("; did not converge in", iter, "iterations")
      }
    }
    cat("\n")
  }
  invisible(x)
}

#' @rdname result-methods
#' @export
print.calibrate_alpha <- function(x, ...) {
  cat(
    "Calibrated intervention on ",
    x$z_name,
    if (!is.null(x$a)) paste0(" with a = ", x$a),
    "\n\n",
    sep = ""
  )
  print(x$calibration)
  cat("\nEstimates under the calibrated intervention:\n")
  print(x$estimate)
  invisible(x)
}

#' @rdname result-methods
#' @export
print.calibration_curve <- function(x, ...) {
  cat(
    "Calibration curve over ",
    length(unique(x$estimate$alpha)),
    " interventions on ",
    x$z_name,
    if (!is.null(x$a)) paste0(" with a = ", x$a),
    "\n\n",
    sep = ""
  )
  print(x$estimate)
  invisible(x)
}

#' @rdname result-methods
#' @export
print.calibrated_contrasts <- function(x, ...) {
  cat(
    "Decomposition of the effect of a on ",
    x$target,
    " into mediated and calibrated contrasts\n\n",
    sep = ""
  )
  print(x$estimate)
  invisible(x)
}

#' @rdname result-methods
#' @export
coef.tmle_alpha <- function(object, ...) {
  stats::setNames(object$estimate$estimate, object$estimate$target)
}

#' @rdname result-methods
#' @export
coef.calibrate_alpha <- coef.tmle_alpha

#' @rdname result-methods
#' @export
confint.tmle_alpha <- function(object, parm, level = 0.95, ...) {
  checkmate::assert_number(level, lower = 0, upper = 1)
  tab <- object$estimate
  if (!missing(parm)) {
    checkmate::assert_subset(parm, tab$target)
    rows <- match(parm, tab$target)
    tab <- tab[rows]
  }
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  out <- tab$estimate + outer(tab$se, stats::qnorm(probs))
  dimnames(out) <- list(
    tab$target,
    paste(format(100 * probs, trim = TRUE, digits = 3), "%")
  )
  out
}

#' @rdname result-methods
#' @export
confint.calibrate_alpha <- confint.tmle_alpha

# "under alpha = 0.5 on z with a = 1", or "without intervention".
.describe_intervention <- function(x) {
  parts <- character()
  if (is.null(x$alpha)) {
    parts <- "under a time-varying intervention (alpha_list)"
  } else if (!is.null(x$z_name) && x$alpha != 1) {
    parts <- paste0("under alpha = ", x$alpha, " on ", x$z_name)
  }
  if (!is.null(x$a)) {
    parts <- c(parts, paste0("with a = ", x$a))
  }
  if (length(parts) == 0) {
    return("without intervention")
  }
  paste(parts, collapse = " ")
}
