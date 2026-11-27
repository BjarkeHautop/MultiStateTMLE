#' Fit the initial hazard models
#'
#' Fits the initial (nuisance) hazard model for each event process (Cox
#' models, or HAL), and the baseline treatment model if requested, and
#' expands the data to the grid of all observed event times up to `tau`. It
#' also builds the discrete product state space of the processes that other
#' processes' hazards depend on. The result is the input to
#' [tmle_estimate()] and [intervention_curve()].
#'
#' @param data A `data.frame` of observed event history data in long format,
#'   with one row per event and columns `id`, `time` and `event` (the name of
#'   the process that jumped), plus any baseline covariates used in
#'   `models`. This is the format of [simevent::simulate_events()].
#' @param models A list with the hazard model of each event process: a
#'   formula such as `death ~ biomarker + illness`, naming the process on the
#'   left and its covariates on the right, or a [hazard_model()] to also restrict
#'   when individuals are at risk of it. Other processes can be used as
#'   covariates by name, giving their number of events so far (see Counts),
#'   and the time of an event through an indicator such as
#'   `I(event_time(illness, 1) > 0.5)` (see [hazard_model()]), also in the
#'   censoring model. The `state_variables` can also be used by name. Every
#'   event in `data` needs a model, including the censoring process.
#' @param tau Positive number. Follow-up horizon; the target parameters are
#'   evaluated at `tau`.
#' @param censoring String. Name of the censoring process in `models`.
#' @param treatment_model Optional logistic regression model for the
#'   baseline treatment, a two-sided formula such as
#'   `treatment ~ biomarker`. Required when `arms` is given.
#' @param arms Integer vector of the baseline treatment values (`0`/`1`) to
#'   prepare interventional hazards and clever weights for. Each value can
#'   then be passed as `a` to [tmle_estimate()] and
#'   [intervention_curve()]. Defaults to both values, `0:1`, when
#'   `treatment_model` is given.
#' @param state_variables A named list of [state_variable()]s: variables
#'   tracked over the event history, such as whether treatment is stopped,
#'   updated by the events of the processes they have rules for. They can be
#'   used by name in `models`, in their `at_risk` conditions and in the
#'   `when` condition of a target in [tmle_estimate()]. Each is a variable
#'   of the state space, with the values it takes in the data up to `tau`,
#'   and the processes it has rules for are in the state space too. Not
#'   supported with `fit = "hal"`.
#' @param control A named list of tuning settings:
#'   * `use_exponential` (`FALSE`): if `TRUE`, always transform the hazard
#'     increments `dH` to probabilities `1 - exp(-dH)`. Otherwise this is
#'     only done if some increments exceed 1.
#' @param hal_control A named list of settings for processes fitted with
#'   `fit = "hal"`, see Details.
#' @param verbose Logical. If `TRUE`, print model fits and progress.
#'
#' @section Parallel computation:
#' The hazards in each state are predicted in parallel on the
#' [mirai::daemons()] set in the session, if any, such as
#' `mirai::daemons(4)` for 4 local processes.
#'
#' @section Counts:
#' A process used in a model or an `at_risk` condition is its number of
#' events so far, and the term is used as written, as in any R model
#' formula. The hazards are computed in a finite state space, so each
#' count is tracked up to the smallest count from which on nothing that
#' uses it changes: 1 for `I(relapse >= 1)`, 3 for `pmin(relapse, 3)` or
#' `factor(pmin(relapse, 3))` (an effect for each count up to 3), and `k`
#' for `event_time(relapse, k)`. A process that can happen at most once, by
#' its `limit` or in the data, is tracked up to 1. A term that keeps
#' changing, such as a plain `relapse` of a recurrent process, a difference
#' of counts such as `I(stop - restart)`, or an `at_risk` condition such as
#' `~ stop == restart`, is an error: use a term that stops changing, or a
#' [state_variable()] for a quantity of several processes. The column
#' `tracked` of `processes` gives the count tracked for each process.
#'
#' The state space holds the states that can be reached from the start
#' under the `at_risk` conditions, and those in the data.
#'
#' @details
#' HAL fits are experimental and require a `fit.hal()` function to be
#' available on the search path. `hal_control` may contain `cut_time` (35),
#' `cut_one_way` (15), `cut_tk` (3), `cut_tk_values` (`NULL`), `two_way`
#' (`NULL`), `reduce_nk` (0.025), `sl` (a list of HAL library
#' specifications, built from `cut_one_way`, `cut_tk`, `cut_tk_values` and
#' `two_way` by default), `screen_two_way` (`FALSE`), `lambda_cvs`,
#' `event_dependent_cv` (`FALSE`), `npenalize_vars` (`NULL`), `v_folds`
#' (10), `seed` (`NULL`), `reduce_seed_dependence` (`FALSE`),
#' `penalize_time` (`FALSE`) and `verbose` (`FALSE`), with defaults in
#' parentheses.
#'
#' Columns of `data` named after a process, such as the event counts that
#' [simevent::simulate_events()] adds, are ignored: the counts are
#' recomputed from `event`.
#'
#' @return An object of class `initial_fit`, a list with elements:
#'   * `n`: the number of individuals.
#'   * `tau`, `treatment_model` and `arms`: as given (`arms` with its
#'     default filled in).
#'   * `models`: the hazard models, as a list of [hazard_model()]s named by
#'     process.
#'   * `processes`: a `data.table` with one row per process and columns
#'     `name`, `type` (`"terminal"`, `"recurrent"` or `"one_jump"`, inferred
#'     from the data), `model`, `at_risk` (whether the process has a
#'     `limit` or `at_risk` condition) and `tracked` (the number of its
#'     events tracked in the state space, see Counts; `NA` if no hazard
#'     depends on it).
#'   * `cox_fits`: the fitted [survival::coxph()] model of each process,
#'     named by process. [as_sim_model()] builds a simulation model from
#'     them.
#'
#'   It also contains the internal elements used by [tmle_estimate()], among
#'   them `tmp_long` (the expanded data with fitted hazards and clever
#'   weights) and `depend_matrix` (the product state space).
#'
#' @examples
#' # Illness-death data: illness is an illness that can occur once, death ends
#' # follow-up, and so does censoring. One row per event:
#' head(illness_death)
#'
#' # One hazard model per process. The death hazard depends on whether the
#' # individual has fallen ill, and the illness can occur at most once:
#' fit <- fit_initial(
#'   illness_death,
#'   models = list(
#'     hazard_model(illness ~ biomarker + treatment, limit = 1),
#'     death ~ biomarker + treatment + illness,
#'     censoring ~ biomarker
#'   ),
#'   tau = 1
#' )
#' fit
#'
#' # The process types are inferred from the data, and the hazard models are
#' # kept for reference:
#' fit$processes
#' fit$models$death
#'
#' # With a baseline treatment, a treatment model prepares estimates under
#' # both arms, for `a` in tmle_estimate():
#' fit_treated <- fit_initial(
#'   illness_death,
#'   models = list(
#'     hazard_model(illness ~ biomarker + treatment, limit = 1),
#'     death ~ biomarker + treatment + illness,
#'     censoring ~ biomarker
#'   ),
#'   tau = 1,
#'   treatment_model = treatment ~ biomarker
#' )
#' fit_treated
#'
#' \donttest{
#' # Data simulated with simevent can be passed as is. Here relapse is
#' # recurrent, and the death hazard grows with the number of relapses so far.
#' # factor(pmin(relapse, 3)) gives an effect of each count up to 3, so 3
#' # relapses are tracked:
#' if (requireNamespace("simevent", quietly = TRUE)) {
#'   model <- simevent::sim_model(
#'     biomarker = simevent::sim_covariate(function(N) rnorm(N)),
#'     relapse = simevent::sim_process("transient", eta = 0.6, nu = 1),
#'     death = simevent::sim_process("terminal", eta = 0.2, nu = 1),
#'     censoring = simevent::sim_process("censoring", eta = 0.2, nu = 1),
#'     effects = list(
#'       simevent::sim_effect(relapse ~ biomarker, coef = 0.5),
#'       simevent::sim_effect(death ~ relapse, coef = 0.4)
#'     )
#'   )
#'   relapses <- simevent::simulate_events(model, n = 200, seed = 1)
#'
#'   fit_initial(
#'     relapses,
#'     models = list(
#'       relapse ~ biomarker,
#'       death ~ biomarker + factor(pmin(relapse, 3)),
#'       censoring ~ biomarker
#'     ),
#'     tau = 2
#'   )
#' }
#' }
#' @export
fit_initial <- function(
  data,
  models,
  tau,
  censoring = "censoring",
  treatment_model = NULL,
  arms = NULL,
  state_variables = list(),
  control = list(),
  hal_control = list(),
  verbose = FALSE
) {
  checkmate::assert_data_frame(data, min.rows = 1)
  checkmate::assert_names(
    names(data),
    must.include = c("id", "time", "event")
  )
  models <- .as_hazard_models(models)
  checkmate::assert_string(censoring)
  if (!censoring %in% names(models)) {
    stop(
      "models has no model for the censoring process \"",
      censoring,
      "\"; give it one, or name the censoring process in `censoring`."
    )
  }

  events <- as.character(data[["event"]])
  if (anyNA(events)) {
    stop("The event column of data has missing values.")
  }
  unmodelled <- setdiff(unique(events), names(models))
  if (length(unmodelled) > 0) {
    stop(
      "data has events without a model in `models`: ",
      paste(unmodelled, collapse = ", "),
      "."
    )
  }

  # Integer codes for the internal `delta` column: 0 for censoring and
  # 1, 2, ... for the other processes, in the order of `models`.
  codes <- stats::setNames(integer(length(models)), names(models))
  others <- names(models) != censoring
  codes[others] <- seq_len(sum(others))

  # Process count columns, as added by simevent, are recomputed from the
  # events.
  data <- as.data.frame(data)
  data <- data[, setdiff(names(data), c(names(models), "delta")), drop = FALSE]
  data[["delta"]] <- unname(codes[events])
  data[["event"]] <- NULL
  .check_state_variables(
    state_variables,
    processes = names(models),
    censoring = censoring,
    columns = names(data)
  )

  # Missing values would be dropped by the Cox fits but carried forward in
  # the expanded data, so refuse them in the columns the models use.
  used <- unique(c(
    "id",
    "time",
    unlist(lapply(models, function(model) all.vars(model$formula[[3]]))),
    if (!is.null(treatment_model)) all.vars(treatment_model)
  ))
  missing <- Filter(
    function(name) anyNA(data[[name]]),
    intersect(used, names(data))
  )
  if (length(missing) > 0) {
    stop(
      "data has missing values in ",
      paste(missing, collapse = ", "),
      ", which the models use; remove or impute them first."
    )
  }

  # event_time() terms become the jump time variables T.<name> of
  # .fit_initial_impl()'s `depend_time`.
  event_times <- .event_time_covariates(
    models,
    times = data[["time"]],
    columns = c(names(data), names(state_variables))
  )
  fit_models <- models
  for (name in names(fit_models)) {
    fit_models[[name]]$formula <- event_times$formulas[[name]]
  }

  fit_types <- Map(
    .hazard_model_to_fit_type,
    fit_models,
    codes,
    MoreArgs = list(state_vars = c(names(models), names(state_variables)))
  )

  fit <- .fit_initial_impl(
    data = data,
    fit_types = fit_types,
    tau = tau,
    treatment_model = treatment_model,
    arms = arms,
    depend_time = event_times$depend_time,
    jump_processes = event_times$jump_processes,
    state_variables = state_variables,
    control = control,
    hal_control = hal_control,
    verbose = verbose
  )
  fit$models <- models
  fit$state_variables <- state_variables
  fit$processes <- data.table(
    name = fit$processes$name,
    type = fit$processes$type,
    model = vapply(
      models,
      function(model) .deparse_one_line(model$formula),
      character(1),
      USE.NAMES = FALSE
    ),
    at_risk = fit$processes$at_risk,
    tracked = fit$processes$tracked
  )
  fit
}

# fit_initial() with the processes given as `fit_types`, a named list of
# list(delta, model, fit, at_risk, limit) for data with an integer `delta`
# column. The `state_variables` are checked by fit_initial().
.fit_initial_impl <- function(
  data,
  fit_types,
  tau,
  treatment_model = NULL,
  arms = NULL,
  depend_time = list(),
  jump_processes = list(),
  state_variables = list(),
  control = list(),
  hal_control = list(),
  verbose = FALSE
) {
  checkmate::assert_data_frame(data, min.rows = 1)
  checkmate::assert_names(
    names(data),
    must.include = c("id", "time", "delta")
  )
  checkmate::assert_list(state_variables, types = "state_variable")
  checkmate::assert_list(jump_processes, types = "character")
  checkmate::assert_subset(names(jump_processes), names(depend_time))
  checkmate::assert_list(fit_types, names = "unique", min.len = 1)
  checkmate::assert_number(tau, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_formula(treatment_model, null.ok = TRUE)
  checkmate::assert_integerish(
    arms,
    lower = 0,
    upper = 1,
    any.missing = FALSE,
    unique = TRUE,
    null.ok = TRUE
  )
  checkmate::assert_list(depend_time, types = "function", names = "unique")
  checkmate::assert_list(control, names = "unique")
  checkmate::assert_subset(names(control), names(.fit_control_defaults))
  control <- utils::modifyList(.fit_control_defaults, control)
  checkmate::assert_flag(
    control$use_exponential,
    .var.name = "control$use_exponential"
  )
  use_exponential <- control$use_exponential
  checkmate::assert_list(hal_control, names = "unique")
  checkmate::assert_subset(names(hal_control), names(.hal_control_defaults))
  checkmate::assert_flag(verbose)

  if (!is.null(treatment_model) && length(treatment_model) != 3) {
    stop("treatment_model must be a two-sided formula, e.g. treatment ~ L0.")
  }
  if (length(arms) > 0 && is.null(treatment_model)) {
    stop("A treatment_model is required when arms is given.")
  }
  if (length(arms) == 0 && !is.null(treatment_model)) {
    arms <- 0:1
  }

  # The inputs as given, to store in the result.
  inputs <- list(
    fit_types = fit_types,
    treatment_model = treatment_model,
    arms = arms
  )

  fit_types <- Map(.normalize_fit_type, fit_types, names(fit_types))

  # The treatment model as a string without whitespace, as the string
  # manipulations on it below expect.
  if (!is.null(treatment_model)) {
    treatment_model <- gsub("\\s+", "", .deparse_one_line(treatment_model))
  }

  hal_control <- utils::modifyList(.hal_control_defaults, hal_control)

  dt <- copy(as.data.table(data))
  n <- uniqueN(dt[["id"]])

  #-- the processes
  process.names <- names(fit_types)
  process.deltas <- vapply(
    fit_types,
    function(fit.type) fit.type[["delta"]],
    numeric(1),
    USE.NAMES = FALSE
  )
  if (anyDuplicated(process.deltas)) {
    stop("Each process in fit_types must have its own delta code.")
  }
  if (sum(process.deltas == 0) != 1) {
    stop("Exactly one process in fit_types must have delta == 0 (censoring).")
  }
  cens.process.id <- which(process.deltas == 0)
  targets <- process.names[-cens.process.id]

  which.hal <- vapply(
    fit_types,
    function(fit.type) fit.type$fit == "hal",
    logical(1)
  )
  any.hal <- any(which.hal)
  if (any.hal && !exists("fit.hal", mode = "function")) {
    stop(
      "fit = \"hal\" requires a fit.hal() function on the search path, ",
      "which MultiStateTMLE does not provide."
    )
  }
  at.risks <- Filter(Negate(is.null), lapply(fit_types, `[[`, "at.risk"))
  process.types <- .process_types(dt, process.deltas, process.names)

  #-- the covariates, and the processes whose counts are covariates (state
  #-- variables); the jump times T.<name> of depend_time are added below.
  jump.specs <- .jump_specs(depend_time, jump_processes)
  T.names <- names(jump.specs)
  # The processes of each jump time are state variables, so the state space
  # knows whether the jump has happened.
  T.processes <- unique(unlist(
    lapply(jump.specs, `[[`, "processes"),
    use.names = FALSE
  ))
  if (any.hal && any(lengths(lapply(jump.specs, `[[`, "processes")) > 1)) {
    stop(
      "event_time() of several processes is not supported with ",
      "fit = \"hal\".",
      call. = FALSE
    )
  }
  if (!all(T.processes %in% targets)) {
    stop(
      "depend_time names jumps of ",
      paste(setdiff(T.processes, targets), collapse = ", "),
      ", which must be processes other than censoring.",
      call. = FALSE
    )
  }
  if (any.hal && length(state_variables) > 0) {
    stop(
      "State variables are not supported with fit = \"hal\".",
      call. = FALSE
    )
  }
  # The processes that state variables have rules for are state variables
  # too, so the state space knows the state after their jumps.
  state.var.names <- names(state_variables)
  rule.processes <- unique(unlist(
    lapply(state_variables, function(variable) names(variable$rules)),
    use.names = FALSE
  ))
  varnames <- union(
    setdiff(.model_variables(fit_types, at.risks), T.names),
    c(T.processes, rule.processes, state.var.names)
  )
  state.processes <- intersect(process.names, varnames)
  state.deltas <- process.deltas[match(state.processes, process.names)]

  A0.name <- NULL
  if (!is.null(treatment_model)) {
    A0.name <- strsplit(treatment_model, "~")[[1]][1]
    varnames <- unique(c(varnames, A0.name))
  }

  #-- the history of each individual at each event
  if (typeof(dt[["time"]]) != "double") {
    warning("The time column is not numeric and will be converted.")
    dt[, time := as.numeric(time)]
  }
  dt[, tstart := c(0, time[-.N]), by = "id"]
  dt[, tstop := time]
  dt <- .add_count_columns(dt, state.processes, state.deltas)
  process.codes <- stats::setNames(as.list(process.deltas), process.names)
  dt <- .add_jump_time_columns(dt, jump.specs, process.codes)
  state.values <- .add_state_variable_columns(
    dt,
    state_variables,
    process.codes,
    tau
  )
  varnames <- c(varnames, T.names)
  state.names <- c(state.processes, state.var.names, T.names)

  #-- the number of events of each process that the state space tracks
  count.upper <- .count_upper(
    state.processes,
    process.types,
    lapply(fit_types, `[[`, "limit")
  )
  needed <- .needed_counts(jump.specs, count.upper)
  count.caps <- .count_caps(
    dt,
    state.processes,
    .count_uses(
      lapply(inputs$fit_types, `[[`, "model"),
      at.risks,
      state.processes
    ),
    count.upper,
    observed = .max_counts(dt, state.processes, state.deltas),
    needed = needed
  )

  #-- the initial fits
  dt[, idN := seq_len(.N), by = "id"]
  if (!is.null(treatment_model)) {
    fit.A0 <- glm(
      as.formula(treatment_model),
      data = dt[idN == 1],
      family = binomial
    )
    if (verbose) .message_print(summary(fit.A0))
  }
  cox <- .fit_cox_models(
    dt,
    fit_types,
    at.risks,
    depend_time,
    verbose
  )

  #-- the data expanded to all event times up to tau
  unique.times <- sort(unique(dt[["time"]]))
  unique.times <- unique.times[unique.times <= tau]
  tmp.long <- .expand_to_grid(
    dt,
    lapply(cox, `[[`, "hazard"),
    varnames,
    unique.times
  )
  tmp.long <- .add_count_columns(tmp.long, state.processes, state.deltas)
  .add_state_variable_columns(tmp.long, state_variables, process.codes, tau)
  tmp.long[, tstart := c(0, time[-.N]), by = "id"]
  tmp.long <- .add_jump_time_columns(
    tmp.long,
    jump.specs,
    process.codes,
    grid = TRUE
  )
  if (length(arms) > 0) {
    tmp.long[["pi.A0.1"]] <- predict(
      fit.A0,
      newdata = tmp.long,
      type = "response"
    )
  }

  hal <- NULL
  if (any.hal) {
    hal <- .hal_fit(
      tmp.long,
      fit_types,
      which.hal,
      at.risks,
      hal_control,
      state.names,
      varnames,
      unique.times,
      A0.name,
      n,
      tau,
      max(c(1L, count.caps)),
      .n_workers(),
      verbose
    )
    tmp.long <- hal$tmp.long
    hal$tmp.long <- NULL
  }

  tmp.long <- tmp.long[time <= tau]
  if (any.hal) {
    hal <- .hal_design(hal, tmp.long, A0.name)
  }

  #-- the observed hazard increments P.<process>, and the censoring survival
  for (name in process.names) {
    if (name %in% names(at.risks)) {
      tmp.long[, (paste0("at.risk.", name)) := at.risks[[name]](tmp.long)]
    }
    P <- paste0("P.", name)
    if (which.hal[[name]]) {
      tmp.long[,
        (P) := as.numeric(.hal_predict(hal, name, tmp.long[["risk.time"]]))
      ]
    } else {
      tmp.long[,
        (P) := .cox_increment(
          cox[[name]][["fit"]],
          tmp.long,
          tmp.long[[paste0("dhazard.", name)]],
          if (name %in% names(at.risks)) tmp.long[[paste0("at.risk.", name)]]
        )
      ]
    }
    if (name == process.names[cens.process.id]) {
      tmp.long[, surv.0 := exp(-cumsum(get(P))), by = "id"]
      tmp.long[, surv.0.1 := c(1, surv.0[-.N]), by = "id"]
    }
  }

  #-- the state space of the state variables
  if (length(state.names) == 0) {
    depend.matrix <- data.table(state = 1L)
  } else {
    depend.matrix <- .state_grid(
      tmp.long,
      state.names,
      count.caps,
      state.values
    )
    if (any.hal) {
      hal.states <- .hal_state_space(hal, tmp.long)
      tmp.long <- hal.states$tmp.long
      if (!is.null(hal.states$depend.matrix)) {
        depend.matrix <- hal.states$depend.matrix
      }
    } else {
      depend.matrix <- .drop_impossible_jump_times(depend.matrix, jump.specs)
      depend.matrix <- .reachable_states(
        depend.matrix,
        observed = .observed_states(tmp.long, state.names, count.caps),
        processes = intersect(state.processes, targets),
        count.caps = count.caps,
        state_variables = state_variables,
        at.risks = at.risks,
        free = T.names
      )
    }
    depend.matrix[, state := seq_len(.N)]
  }

  #-- clever weights, and the hazards in each state under each arm
  intervene.A0 <- length(arms) > 0
  if (intervene.A0) {
    setnames(tmp.long, A0.name, "A0.obs")
  }
  tmp.long[,
    C := cumsum(1 * (time == time.obs & delta == 0 & time %in% dt[["time"]])),
    by = "id"
  ]
  tmp.long[, C.1 := c(0, C[-.N]), by = "id"]
  tmp.long[, clever.weight := (C.1 == 0) / surv.0.1]
  if (any.hal) {
    tmp.long[, state := 0]
  } else {
    tmp.long[,
      state := .match_states(tmp.long, depend.matrix, state.names, count.caps)
    ]
  }

  for (aa in if (intervene.A0) arms else 1) {
    tmp.long <- .predict_states(
      tmp.long,
      aa,
      intervene.A0,
      A0.name,
      targets,
      cox,
      hal,
      at.risks,
      depend.matrix,
      varnames,
      verbose
    )
  }

  for (name.P in names(tmp.long)[startsWith(names(tmp.long), "P.")]) {
    if (use_exponential || any(tmp.long[[name.P]] > 1)) {
      if (!use_exponential && verbose) {
        message(
          "Transforming ",
          name.P,
          " with 1 - exp(-x) to avoid values > 1."
        )
      }
      tmp.long[, (name.P) := 1 - exp(-tmp.long[[name.P]])]
    }
  }

  #-- the jump time rows of the state space, and the links between states
  var.depend <- setdiff(names(depend.matrix), "state")
  process.depend <- intersect(var.depend, process.names)
  if (any.hal) {
    link.specs <- .jump_specs_from_names(.hal_time_vars(hal))
  } else {
    # The jump times of the processes in the state space, ordered by their
    # first process.
    first <- vapply(jump.specs, function(spec) spec$processes[1], character(1))
    in.states <- vapply(
      jump.specs,
      function(spec) all(spec$processes %in% process.depend),
      logical(1)
    )
    order.specs <- order(match(first, process.depend))
    link.specs <- jump.specs[order.specs[in.states[order.specs]]]
  }
  depend.matrix <- .add_jump_time_rows(
    depend.matrix,
    tmp.long,
    link.specs,
    zero = any.hal
  )
  linked <- .link_states(
    depend.matrix,
    tmp.long,
    var.depend,
    process.depend,
    link.specs,
    state_variables,
    if (any.hal) NULL else count.caps
  )

  processes <- data.table(
    name = process.names,
    delta = process.deltas,
    type = unlist(process.types, use.names = FALSE),
    model = vapply(
      inputs$fit_types,
      function(spec) .deparse_one_line(spec$model),
      character(1),
      USE.NAMES = FALSE
    ),
    at_risk = process.names %in% names(at.risks),
    tracked = unname(count.caps[process.names])
  )

  structure(
    list(
      n = n,
      tau = tau,
      processes = processes,
      treatment_model = inputs$treatment_model,
      arms = inputs$arms,
      tmp_long = linked$tmp.long,
      depend_matrix = linked$depend.matrix,
      process_names = process.names,
      process_deltas = process.deltas,
      process_types = process.types,
      cens_process_id = cens.process.id,
      at_risks = at.risks,
      count_caps = count.caps,
      count_upper = count.upper,
      cox_fits = lapply(cox, function(x) .detach_cox_env(x$fit))
    ),
    class = "initial_fit"
  )
}
