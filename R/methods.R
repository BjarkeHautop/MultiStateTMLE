#' Methods for the result objects
#'
#' Print, extract the estimates of, and compute confidence intervals for the
#' objects returned by [fit_initial()], [tmle_estimate()],
#' [intervention_curve()] and [calibrated_contrasts()].
#'
#' @param x,object A result object.
#' @param parm Optional character vector of target processes to compute
#'   confidence intervals for. Defaults to all.
#' @param level Number between 0 and 1. Confidence level.
#' @param ... Ignored.
#'
#' @return `print()` returns `x` invisibly. `coef()` returns a numeric vector
#'   of estimates named by target process, and `confint()` a matrix of Wald
#'   confidence intervals with one row per target process.
#'
#' @examples
#' fit <- fit_initial(
#'   illness_death,
#'   models = list(
#'     hazard_model(z ~ L0, limit = 1),
#'     death ~ L0 + z,
#'     censoring ~ L0
#'   ),
#'   tau = 1
#' )
#' fit
#'
#' est <- tmle_estimate(fit, risk(death, z), scale_intensity(z = 0.5))
#' est
#'
#' # The estimates, named by target process:
#' coef(est)
#'
#' # Confidence intervals at another level, for all or some targets:
#' confint(est, level = 0.9)
#' confint(est, "death")
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
print.tmle_estimate <- function(x, ...) {
  method <- if (x$method == "one_step") "One-step" else "TMLE"
  cat(
    method,
    " estimate of ",
    format(x$target),
    "
",
    sep = ""
  )
  cat(
    .describe_intervention(x),
    "

",
    sep = ""
  )
  if (.is_calibrated(x)) {
    print(x$calibration)
    cat(
      "
Estimates under the calibrated intervention:
"
    )
    print(x$estimate)
    return(invisible(x))
  }
  print(x$estimate)
  cat(
    "
"
  )
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
    cat(
      "
"
    )
  }
  invisible(x)
}

#' @rdname result-methods
#' @export
print.intervention_curve <- function(x, ...) {
  cat(
    "Intervention curve of ",
    format(x$target),
    " against ",
    format(x$exposure),
    "
over ",
    length(unique(x$estimate$alpha)),
    " scalings of the intensity of ",
    x$exposure$processes,
    if (!is.null(x$a)) paste0(", with a = ", x$a),
    "

",
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
    " into mediated and calibrated contrasts

",
    sep = ""
  )
  print(x$estimate)
  invisible(x)
}

#' @rdname result-methods
#' @export
coef.tmle_estimate <- function(object, ...) {
  stats::setNames(object$estimate$estimate, object$estimate$target)
}

#' @rdname result-methods
#' @export
confint.tmle_estimate <- function(object, parm, level = 0.95, ...) {
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

# "under scale_intensity(z = 0.5), with a = 1", or "without intervention".
.describe_intervention <- function(x) {
  parts <- character()
  if (!is.null(x$intervention)) {
    parts <- paste("under", format(x$intervention))
  }
  if (!is.null(x$a)) {
    parts <- c(parts, paste0("with a = ", x$a))
  }
  if (length(parts) == 0) {
    return("without intervention")
  }
  paste(parts, collapse = ", ")
}
