#' Define a state variable
#'
#' Describes a variable that [fit_initial()] tracks over the event history
#' of each individual, for use in the hazard models. It starts at `start`,
#' and each event of a process named in `...` updates it by that process's
#' rule. Pass it to [fit_initial()] in `state_variables`, named by the
#' variable, which can then be used by name in the `models` formulas and
#' their `at_risk` conditions, and in the `when` condition of a target in
#' [tmle_estimate()].
#'
#' A state variable is a state of its own, tracked exactly. Use one for a
#' quantity of several processes, such as whether treatment is stopped
#' (stops minus restarts): in a model it would keep changing with the
#' counts, which the state space of [fit_initial()] only tracks up to the
#' count from which on the models stop changing (see its Counts section).
#'
#' @param start Number. The value of the variable before any event.
#' @param ... The rules of the variable, named by the process whose events
#'   they apply to: a number, the value after an event of the process, or a
#'   vectorised function of the value before the event that returns the
#'   value after it, such as `\(s) pmin(s + 1, 3)`. Events of other
#'   processes leave the variable unchanged.
#'
#' @return An object of class `state_variable`.
#'
#' @examples
#' # Whether treatment is stopped: stop events (stop) set it to 1, and
#' # restart events (restart) set it back to 0:
#' state_variable(start = 0, stop = 1, restart = 0)
#'
#' # The number of stops, counted up to 3:
#' state_variable(start = 0, stop = \(s) pmin(s + 1, 3))
#'
#' \donttest{
#' # Treatment that is stopped and restarted, any number of times, with a
#' # death hazard that is higher while off treatment:
#' if (requireNamespace("simevent", quietly = TRUE)) {
#'   model <- simevent::sim_model(
#'     biomarker = simevent::sim_covariate(function(N) rnorm(N)),
#'     stop = simevent::sim_process("transient", eta = 0.8, nu = 1),
#'     restart = simevent::sim_process("transient", eta = 0.8, nu = 1),
#'     death = simevent::sim_process("terminal", eta = 0.2, nu = 1),
#'     censoring = simevent::sim_process("censoring", eta = 0.2, nu = 1),
#'     effects = list(
#'       simevent::sim_effect(stop ~ biomarker, coef = 0.5),
#'       simevent::sim_effect(death ~ (stop - restart), coef = 1),
#'       # Stops only while on treatment, restarts only while off it:
#'       simevent::sim_effect(stop ~ (stop != restart), coef = -1e5),
#'       simevent::sim_effect(restart ~ (restart != stop - 1), coef = -1e5)
#'     )
#'   )
#'   treatment <- simevent::simulate_events(model, n = 200, seed = 1)
#'
#'   # Some individuals stop several times. Off treatment as I(stop - restart)
#'   # keeps changing with the counts, which the state space can't track,
#'   # while the state variable off is exact:
#'   fit <- fit_initial(
#'     treatment,
#'     models = list(
#'       hazard_model(stop ~ biomarker, at_risk = ~ off == 0),
#'       hazard_model(restart ~ biomarker, at_risk = ~ off == 1),
#'       death ~ biomarker + off,
#'       censoring ~ biomarker
#'     ),
#'     tau = 2,
#'     state_variables = list(
#'       off = state_variable(start = 0, stop = 1, restart = 0)
#'     )
#'   )
#'   print(fit)
#'
#'   # A stop leads to a state with off = 1, and a restart to one with
#'   # off = 0:
#'   print(fit$depend_matrix)
#'
#'   tmle_estimate(fit, risk(death))
#' }
#' }
#' @export
state_variable <- function(start, ...) {
  checkmate::assert_number(start, finite = TRUE)
  rules <- list(...)
  if (length(rules) == 0) {
    stop(
      "A state variable needs a rule for at least one process, ",
      "e.g. state_variable(start = 0, stop = 1).",
      call. = FALSE
    )
  }
  if (is.null(names(rules)) || any(names(rules) == "")) {
    stop(
      "Each rule of a state variable must be named by its process, ",
      "e.g. state_variable(start = 0, stop = 1).",
      call. = FALSE
    )
  }
  if (anyDuplicated(names(rules))) {
    stop(
      "A state variable can have one rule per process, but ",
      paste(unique(names(rules)[duplicated(names(rules))]), collapse = ", "),
      " has several.",
      call. = FALSE
    )
  }
  for (process in names(rules)) {
    rule <- rules[[process]]
    if (
      !is.function(rule) &&
        !checkmate::test_number(rule, finite = TRUE)
    ) {
      stop(
        "The rule for ",
        process,
        " must be a number or a function of the value before the event.",
        call. = FALSE
      )
    }
  }
  structure(
    list(start = as.numeric(start), rules = rules),
    class = "state_variable"
  )
}

#' @export
print.state_variable <- function(x, ...) {
  cat("State variable starting at ", format(x$start), "\n", sep = "")
  for (process in names(x$rules)) {
    rule <- x$rules[[process]]
    cat(
      "  after ",
      process,
      ": ",
      if (is.function(rule)) .deparse_one_line(rule) else format(rule),
      "\n",
      sep = ""
    )
  }
  invisible(x)
}

# Column names used internally by fit_initial() and tmle_estimate(), and
# their prefixes, that a state variable may not take.
.reserved_names <- c(
  "id",
  "time",
  "event",
  "delta",
  "state",
  "tstart",
  "tstop",
  "idN",
  "row_id",
  "time.obs",
  "final.time",
  "grid.time",
  "risk.time",
  "surv.0",
  "surv.0.1",
  "C",
  "C.1",
  "Q",
  "at.risk",
  "state.row.index",
  "pi.A0.1",
  "A0.obs",
  "target.only.in.state",
  "clever.weight"
)
.reserved_prefixes <- c(
  "T.",
  "time.",
  "s.",
  "P.",
  "Q.",
  "at.risk.",
  "dhazard.",
  "hazard.",
  "clever.",
  "alpha."
)

# Check fit_initial()'s `state_variables` against the names of the
# processes and the data columns.
.check_state_variables <- function(
  state_variables,
  processes,
  censoring,
  columns
) {
  checkmate::assert_list(
    state_variables,
    types = "state_variable",
    names = "unique"
  )
  for (name in names(state_variables)) {
    what <- paste0("The state variable ", name)
    if (!identical(make.names(name), name)) {
      stop(what, " must have a syntactic name.", call. = FALSE)
    }
    if (name %in% processes) {
      stop(what, " has the name of a process.", call. = FALSE)
    }
    if (name %in% columns) {
      stop(what, " has the name of a column of data.", call. = FALSE)
    }
    if (
      name %in% .reserved_names || any(startsWith(name, .reserved_prefixes))
    ) {
      stop(what, " has a name used internally; rename it.", call. = FALSE)
    }
    rule_processes <- names(state_variables[[name]]$rules)
    unknown <- setdiff(rule_processes, processes)
    if (length(unknown) > 0) {
      stop(
        what,
        " has rules for ",
        paste(unknown, collapse = ", "),
        ", which have no model in `models`.",
        call. = FALSE
      )
    }
    if (censoring %in% rule_processes) {
      stop(
        what,
        " has a rule for the censoring process ",
        censoring,
        ", which ends follow-up.",
        call. = FALSE
      )
    }
  }
  invisible(state_variables)
}

# The value of state variable `variable` (named `name`) after an event of
# `process`, from the values `values` before it.
.apply_rule <- function(variable, name, process, values) {
  rule <- variable$rules[[process]]
  if (!is.function(rule)) {
    return(rep(rule, length(values)))
  }
  after <- rule(values)
  if (
    !is.numeric(after) ||
      length(after) != length(values) ||
      !all(is.finite(after))
  ) {
    stop(
      "The rule for ",
      process,
      " of the state variable ",
      name,
      " must return a finite number for each value, as a vector of the ",
      "same length.",
      call. = FALSE
    )
  }
  as.numeric(after)
}

# Add to `x` (sorted by id and time, with `delta`) the value of each state
# variable of `state_variables` before the event of each row, replaying the
# events of each individual through the rules. `codes` are the delta codes
# of the processes, named by process. Adds the columns by reference, and
# returns the values each variable takes up to time `tau`, named by
# variable.
.add_state_variable_columns <- function(x, state_variables, codes, tau) {
  values <- list()
  n <- nrow(x)
  first <- !duplicated(x[["id"]])
  id.index <- cumsum(first)
  for (name in names(state_variables)) {
    variable <- state_variables[[name]]
    rule.codes <- unlist(codes[names(variable$rules)], use.names = FALSE)
    rule.processes <- names(variable$rules)
    events <- which(x[["delta"]] %in% rule.codes)
    # The value after each event, from the value after the previous event
    # of the same individual.
    after <- rep(NA_real_, n)
    current <- rep(variable$start, max(c(0, id.index)))
    order.in.id <- stats::ave(events, id.index[events], FUN = seq_along)
    for (k in seq_len(max(c(0, order.in.id)))) {
      rows <- events[order.in.id == k]
      value <- current[id.index[rows]]
      delta <- x[["delta"]][rows]
      for (j in seq_along(rule.processes)) {
        jumps <- delta == rule.codes[j]
        if (any(jumps)) {
          value[jumps] <- .apply_rule(
            variable,
            name,
            rule.processes[j],
            value[jumps]
          )
        }
      }
      current[id.index[rows]] <- value
      after[rows] <- value
    }
    # The value before each row is the value after the last event before it.
    before <- c(NA_real_, after[-n])
    before[first] <- variable$start
    before <- .locf(before)
    set(x, j = name, value = before)
    values[[name]] <- sort(unique(c(
      variable$start,
      after[events[x[["time"]][events] <= tau]]
    )))
  }
  invisible(values)
}
