#' Target parameters
#'
#' Describe what [tmle_estimate()] and [intervention_curve()] estimate for
#' one or more event processes, by the horizon `tau` of the initial fit:
#' * `risk()`: the risk of an event by `tau`, \eqn{E[N(\tau)]}, which is the
#'   expected number of events for a recurrent process.
#' * `years_lost()`: the expected time lost before `tau`, the integral of
#'   the risk up to `tau`.
#'
#' @param ... The event processes, as bare names (`death`) or strings
#'   (`"death"`). Each is estimated separately, under the same intervention.
#'   To give the names in a character vector `x`, use
#'   `do.call(risk, as.list(x))`.
#' @param when Optional one-sided formula of a condition on the state, such
#'   as `~ z == 0`: only events while it holds count. It may use the
#'   processes that hazards depend on (their number of events so far) and
#'   the `derived_vars` of [fit_initial()].
#' @param only_first Logical. If `TRUE`, only the first event of each
#'   process counts.
#'
#' @return An object of class `tmle_target`.
#'
#' @examples
#' # The risk of death by tau, and the expected years lost before it:
#' risk(death)
#' years_lost(death)
#'
#' # Several processes, each estimated separately under the same intervention:
#' risk(death, z)
#'
#' # Only deaths before illness count:
#' risk(death, when = ~ z == 0)
#'
#' # For a recurrent process, the risk of a first event rather than the
#' # expected number of events:
#' risk(relapse, only_first = TRUE)
#'
#' # Processes held in a character vector:
#' targets <- c("death", "z")
#' do.call(risk, as.list(targets))
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
#' tmle_estimate(fit, risk(death, z))
#' tmle_estimate(fit, risk(death, when = ~ z == 0))$estimate
#' @name target-parameters
NULL

#' @rdname target-parameters
#' @export
risk <- function(..., when = NULL, only_first = FALSE) {
  .tmle_target(
    .process_names(substitute(list(...)), "risk"),
    "risk",
    when,
    only_first
  )
}

#' @rdname target-parameters
#' @export
years_lost <- function(..., when = NULL, only_first = FALSE) {
  .tmle_target(
    .process_names(substitute(list(...)), "years_lost"),
    "years_lost",
    when,
    only_first
  )
}

.tmle_target <- function(processes, estimand, when, only_first) {
  checkmate::assert_formula(when, null.ok = TRUE)
  if (!is.null(when) && length(when) != 2) {
    stop("`when` must be a one-sided formula, e.g. ~ z == 0.")
  }
  checkmate::assert_flag(only_first)
  structure(
    list(
      processes = processes,
      estimand = estimand,
      when = when,
      only_first = only_first
    ),
    class = "tmle_target"
  )
}

# The process names in the unevaluated `list(...)` call `args`, each a bare
# name or a string.
.process_names <- function(args, fun) {
  args <- as.list(args)[-1]
  if (length(args) == 0) {
    stop(fun, "() needs at least one process.")
  }
  if (!is.null(names(args)) && any(names(args) != "")) {
    stop(
      "Unknown argument(s) to ",
      fun,
      "(): ",
      paste(names(args)[names(args) != ""], collapse = ", "),
      "."
    )
  }
  processes <- vapply(
    args,
    function(arg) {
      if (is.name(arg)) {
        return(as.character(arg))
      }
      if (is.character(arg) && length(arg) == 1) {
        return(arg)
      }
      stop(
        "The processes in ",
        fun,
        "() must be bare names or strings, not ",
        .deparse_one_line(arg),
        "."
      )
    },
    character(1)
  )
  if (anyDuplicated(processes)) {
    stop("Each process can only be given once to ", fun, "().")
  }
  processes
}

#' @export
format.tmle_target <- function(x, ...) {
  args <- x$processes
  if (!is.null(x$when)) {
    args <- c(args, paste("when =", .deparse_one_line(x$when)))
  }
  if (x$only_first) {
    args <- c(args, "only_first = TRUE")
  }
  paste0(x$estimand, "(", paste(args, collapse = ", "), ")")
}

#' @export
print.tmle_target <- function(x, ...) {
  cat("Target: ", format(x), "\n", sep = "")
  invisible(x)
}

# `target` as a tmle_target with processes in `initial_fit`.
.assert_target <- function(target, initial_fit, name = "target") {
  if (!inherits(target, "tmle_target")) {
    stop("`", name, "` must be given by risk() or years_lost().")
  }
  unknown <- setdiff(target$processes, initial_fit$process_names)
  if (length(unknown) > 0) {
    stop(
      "The process(es) ",
      paste(unknown, collapse = ", "),
      " in `",
      name,
      "` are not in initial_fit, which has: ",
      paste(initial_fit$process_names, collapse = ", "),
      "."
    )
  }
  invisible(target)
}
