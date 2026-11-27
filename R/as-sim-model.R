#' Build a simulation model from an initial fit
#'
#' Turns the Cox models of an [fit_initial()] fit into a
#' [simevent::sim_model()], with [simevent::sim_model_from_fits()], so data
#' resembling the observed data can be simulated with
#' [simevent::simulate_events()], e.g. to check an estimator where the truth
#' is known. The process types and limits are taken from the fit.
#'
#' Covariates are resampled from `data`, one row per individual. A process
#' that a hazard depends on enters it as in the fit: a recurrent process as
#' the indicators of at least 1, ..., `max_count` jumps, and any other as its
#' number of jumps so far. A process that jumps at most once in the data
#' (type `"one_jump"`) can jump at most once in the simulation too.
#'
#' Fits with an `at_risk` condition in a [hazard_model()], with
#' `event_time()` covariates, or with `state_variables`, are not supported.
#'
#' @param fit An initial fit, as returned by [fit_initial()].
#' @param data The `data` that `fit` was fitted to, to resample covariates
#'   from.
#'
#' @return A [simevent::sim_model()].
#' @examples
#' \donttest{
#' if (requireNamespace("simevent", quietly = TRUE)) {
#'   fit <- fit_initial(
#'     illness_death,
#'     models = list(
#'       hazard_model(illness ~ biomarker + treatment, limit = 1),
#'       death ~ biomarker + treatment + illness,
#'       censoring ~ biomarker
#'     ),
#'     tau = 1
#'   )
#'   model <- as_sim_model(fit, illness_death)
#'   simulated <- simevent::simulate_events(model, n = 500, seed = 1)
#'   head(simulated)
#' }
#' }
#' @export
as_sim_model <- function(fit, data) {
  checkmate::assert_class(fit, "initial_fit")
  checkmate::assert_data_frame(data, min.rows = 1)
  if (!requireNamespace("simevent", quietly = TRUE)) {
    stop("as_sim_model() requires the simevent package.")
  }

  processes <- fit$processes
  custom <- processes$name[vapply(
    fit$models[processes$name],
    function(model) !is.null(model$at_risk),
    logical(1)
  )]
  if (length(custom) > 0) {
    stop(
      "as_sim_model() can't simulate a process with an `at_risk` ",
      "condition, which simevent has no counterpart for: ",
      paste(custom, collapse = ", "),
      "."
    )
  }
  if (length(fit$state_variables) > 0) {
    stop(
      "as_sim_model() can't simulate a fit with state variables (",
      paste(names(fit$state_variables), collapse = ", "),
      "), which simevent has no counterpart for."
    )
  }
  # Variables of the fits that are neither data columns nor processes are
  # the jump time variables T.<name> of event_time().
  fit_vars <- unique(unlist(lapply(
    fit$cox_fits,
    function(cox) all.vars(stats::delete.response(stats::terms(cox)))
  )))
  unknown <- setdiff(fit_vars, c(names(data), processes$name))
  if (length(unknown) > 0) {
    stop(
      "as_sim_model() can't simulate a fit with event_time() covariates (",
      paste(unknown, collapse = ", "),
      "), which simevent has no counterpart for."
    )
  }

  censoring <- processes$name[fit$process_deltas == 0]
  types <- stats::setNames(
    ifelse(processes$type == "terminal", "terminal", "transient"),
    processes$name
  )
  types[censoring] <- "censoring"
  limits <- unlist(lapply(
    stats::setNames(nm = processes$name[types == "transient"]),
    function(name) {
      limit <- fit$models[[name]]$limit
      if (
        is.null(limit) && processes$type[processes$name == name] == "one_jump"
      ) {
        limit <- 1
      }
      limit
    }
  ))
  if (is.null(limits)) {
    limits <- numeric()
  }

  simevent::sim_model_from_fits(fit$cox_fits, data, types, limits)
}
