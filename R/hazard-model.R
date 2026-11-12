#' Specify the hazard model of an event process
#'
#' Describes how [fit_initial()] models the hazard of one event process. A
#' plain formula such as `death ~ L0 + z` in `models` is short for
#' `hazard_model(death ~ L0 + z)`; use `hazard_model()` to also restrict
#' when an individual is at risk of the process, or to change how it is
#' fitted.
#'
#' @param formula A two-sided formula. The left-hand side names the process,
#'   matching its name in the `event` column of the data, and the right-hand
#'   side gives the covariates of its hazard. Other processes can be used as
#'   covariates by name, giving their number of events so far.
#' @param limit Optional positive integer. Maximum number of events of the
#'   process per individual: individuals are only at risk while they have
#'   had fewer than `limit` events. `limit = 1` is an event that can happen
#'   at most once, as `sim_process(limit = 1)` in simevent.
#' @param at_risk Optional one-sided formula, such as `~ z == 0`, of a
#'   condition on the processes (their number of events so far) and the
#'   `derived_vars` of [fit_initial()]. Individuals are only at risk of the
#'   process while it holds. Combined with `limit` if both are given.
#' @param fit How the hazard is fitted: `"cox"` (default), a Cox model, or
#'   `"hal"`, the highly adaptive lasso (experimental, see [fit_initial()]).
#'
#' @return An object of class `hazard_model`.
#'
#' @examples
#' # A plain formula in fit_initial()'s `models` is short for
#' # hazard_model(formula):
#' hazard_model(death ~ L0 + z)
#'
#' # An illness that can happen at most once:
#' hazard_model(z ~ L0, limit = 1)
#'
#' # A relapse that can happen at most twice, and only before death of
#' # another cause, w:
#' hazard_model(relapse ~ L0, limit = 2, at_risk = ~ w == 0)
#'
#' # Use them in fit_initial(), mixed with plain formulas:
#' fit <- fit_initial(
#'   illness_death,
#'   models = list(
#'     hazard_model(z ~ L0, limit = 1),
#'     death ~ L0 + z,
#'     censoring ~ L0
#'   ),
#'   tau = 1
#' )
#' fit$processes
#' @export
hazard_model <- function(
  formula,
  limit = NULL,
  at_risk = NULL,
  fit = c("cox", "hal")
) {
  checkmate::assert_formula(formula)
  if (length(formula) != 3 || !is.name(formula[[2]])) {
    stop(
      "`formula` must be two-sided with the process name on the left, ",
      "e.g. death ~ L0 + z."
    )
  }
  checkmate::assert_count(limit, positive = TRUE, null.ok = TRUE)
  checkmate::assert_formula(at_risk, null.ok = TRUE)
  if (!is.null(at_risk) && length(at_risk) != 2) {
    stop("`at_risk` must be a one-sided formula, e.g. ~ z == 0.")
  }
  fit <- match.arg(fit)
  structure(
    list(
      name = as.character(formula[[2]]),
      formula = formula,
      limit = limit,
      at_risk = at_risk,
      fit = fit
    ),
    class = "hazard_model"
  )
}

#' @export
print.hazard_model <- function(x, ...) {
  cat("Hazard model: ", .deparse_one_line(x$formula), "\n", sep = "")
  if (!is.null(x$limit)) {
    cat("  limit: ", x$limit, "\n", sep = "")
  }
  if (!is.null(x$at_risk)) {
    cat("  at risk while: ", .deparse_one_line(x$at_risk[[2]]), "\n", sep = "")
  }
  if (x$fit != "cox") {
    cat("  fit: ", x$fit, "\n", sep = "")
  }
  invisible(x)
}

# fit_initial()'s `models` as a list of hazard_model objects named by
# process.
.as_hazard_models <- function(models) {
  if (inherits(models, c("formula", "hazard_model"))) {
    models <- list(models)
  }
  checkmate::assert_list(models, min.len = 1)
  given_names <- names(models)
  models <- lapply(seq_along(models), function(i) {
    model <- models[[i]]
    if (inherits(model, "formula")) {
      model <- hazard_model(model)
    }
    if (!inherits(model, "hazard_model")) {
      stop("models[[", i, "]] must be a formula or a hazard_model().")
    }
    if (!is.null(given_names) && !given_names[i] %in% c("", model$name)) {
      stop(
        "models[[",
        i,
        "]] is named \"",
        given_names[i],
        "\" but models the process \"",
        model$name,
        "\"."
      )
    }
    model
  })
  names(models) <- vapply(models, function(model) model$name, character(1))
  duplicated_names <- unique(names(models)[duplicated(names(models))])
  if (length(duplicated_names) > 0) {
    stop(
      "Each process must have one model, but ",
      paste(duplicated_names, collapse = ", "),
      " has several."
    )
  }
  models
}

# A hazard_model as an entry of the internal `fit_types`,
# list(delta, model, fit, at_risk), where `at_risk` is a function of the
# data as .fit_initial_impl() expects. `state_vars` are the processes and
# derived variables that `at_risk` may use.
.hazard_model_to_fit_type <- function(model, delta, state_vars) {
  env <- environment(model$formula)
  out <- list(
    delta = delta,
    model = stats::as.formula(call("~", model$formula[[3]]), env = env),
    fit = model$fit
  )
  conditions <- list()
  if (!is.null(model$limit)) {
    conditions <- c(
      conditions,
      list(call("<", as.name(model$name), model$limit))
    )
  }
  if (!is.null(model$at_risk)) {
    env <- environment(model$at_risk)
    conditions <- c(conditions, list(model$at_risk[[2]]))
  }
  if (length(conditions) > 0) {
    out$at_risk <- .state_function(
      Reduce(function(x, y) call("&", x, y), conditions),
      state_vars,
      env,
      paste("The at-risk condition of", model$name)
    )
  }
  out
}

# A function of a data.table `dt` (the data or the state space) evaluating
# `expr`, with the state variables in `expr` read as `dt[["<name>"]]` and
# other variables from `env`. Writing the state variables out this way lets
# .fit_initial_impl() find the ones an `at_risk` function uses. `what`
# describes `expr` in error messages.
.state_function <- function(expr, state_vars, env, what) {
  vars <- all.vars(expr)
  unknown <- vars[
    !vars %in% state_vars &
      !vapply(vars, exists, logical(1), envir = env)
  ]
  if (length(unknown) > 0) {
    stop(
      what,
      " uses ",
      paste(unknown, collapse = ", "),
      ", but it may only use the processes and derived_vars that hazards ",
      "depend on.",
      call. = FALSE
    )
  }
  in_state <- vars[vars %in% state_vars]
  lookups <- lapply(in_state, function(var) call("[[", quote(dt), var))
  names(lookups) <- in_state
  fun <- function(dt) NULL
  body(fun) <- do.call(substitute, list(expr, lookups))
  environment(fun) <- env
  fun
}
