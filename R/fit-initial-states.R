# The discrete state space of fit_initial(): the product of the values of the
# state variables, the hazards predicted in each state, and the links from
# each state to the state after a jump.

# All combinations of the values of the state variables: the counts 0 to
# `count.caps` of each process, and the values in `values` of each of
# fit_initial()'s `state_variables`. A jump time variable T.<name> takes the
# values it has in the data and, from `time.<name>`, the values it would get
# from a jump at any time of the grid.
.state_grid <- function(tmp.long, state.names, count.caps, values = list()) {
  grid.list <- lapply(state.names, function(varname) {
    if (varname %in% names(values)) {
      return(values[[varname]])
    }
    if (varname %in% names(count.caps)) {
      return(0:count.caps[[varname]])
    }
    time.col <- sub("^T\\.", "time.", varname)
    sort(unique(c(tmp.long[[varname]], tmp.long[[time.col]])))
  })
  depend.matrix <- as.data.table(do.call(expand.grid, grid.list))
  setnames(depend.matrix, state.names)
  depend.matrix
}

# Drop the states in which the jump time variable of the kth jump of its
# processes (`specs`, see .jump_specs()) differs from its value before the
# jump, its function at 0, while the processes have made fewer than k jumps.
.drop_impossible_jump_times <- function(depend.matrix, specs) {
  for (T.var in names(specs)) {
    spec <- specs[[T.var]]
    before <- as.numeric(spec$fun(0))
    depend.matrix <- depend.matrix[
      !(.jump_count(depend.matrix, spec) <= spec$k - 1 &
        depend.matrix[[T.var]] != before)
    ]
  }
  depend.matrix
}

# The state variables of each row of `x` as they are in the state space:
# the counts capped at `count.caps`, the other variables as they are.
.capped_states <- function(x, state.names, count.caps) {
  as.data.table(lapply(
    stats::setNames(nm = state.names),
    function(name) {
      value <- as.numeric(x[[name]])
      if (name %in% names(count.caps)) {
        value <- pmin(value, count.caps[[name]])
      }
      value
    }
  ))
}

# The distinct states of the rows of tmp.long.
.observed_states <- function(tmp.long, state.names, count.caps) {
  unique(.capped_states(tmp.long, state.names, count.caps))
}

# The state of each row of tmp.long, matched exactly by its state
# variables in the state space `depend.matrix` (one row per state).
.match_states <- function(tmp.long, depend.matrix, state.names, count.caps) {
  if (length(state.names) == 0) {
    return(rep(1L, nrow(tmp.long)))
  }
  observed <- .capped_states(tmp.long, state.names, count.caps)
  states <- .capped_states(depend.matrix, state.names, count.caps)
  states[, state := depend.matrix[["state"]]]
  matched <- states[observed, on = state.names, state, mult = "first"]
  if (anyNA(matched)) {
    stop(
      "Internal error: rows of the data in no state of the state space.",
      call. = FALSE
    )
  }
  matched
}

# Keep the states of `depend.matrix` that can be reached from the start (no
# events, and each state variable at its start), or that are `observed`.
# From a state, each of `processes` at risk in it (by `at.risks`) leads to
# the state with its count one more (capped at `count.caps`) and the state
# variables with a rule for it set by the rule. The variables in `free` (the
# jump time variables) are left out: all their values are kept for the
# states reached.
.reachable_states <- function(
  depend.matrix,
  observed,
  processes,
  count.caps,
  state_variables,
  at.risks,
  free = character()
) {
  key.names <- setdiff(names(depend.matrix), c("state", free))
  if (length(key.names) == 0) {
    return(depend.matrix)
  }
  as.key <- function(x) {
    columns <- lapply(key.names, function(name) as.numeric(x[[name]]))
    do.call(paste, c(columns, sep = "|"))
  }
  grid <- unique(depend.matrix[, key.names, with = FALSE])
  grid.keys <- as.key(grid)

  start <- as.list(stats::setNames(numeric(length(key.names)), key.names))
  for (name in intersect(key.names, names(state_variables))) {
    start[[name]] <- state_variables[[name]]$start
  }
  reached <- unique(c(as.key(start), as.key(observed)))
  frontier <- grid[grid.keys %in% reached]
  while (nrow(frontier) > 0) {
    next.keys <- character()
    for (process in processes) {
      at.risk <- rep(TRUE, nrow(frontier))
      if (process %in% names(at.risks)) {
        at.risk <- rep_len(
          as.logical(at.risks[[process]](frontier)),
          nrow(frontier)
        )
      }
      if (!any(at.risk)) {
        next
      }
      before <- frontier[at.risk]
      after <- copy(before)
      if (process %in% key.names) {
        set(
          after,
          j = process,
          value = pmin(after[[process]] + 1, count.caps[[process]])
        )
      }
      for (name in names(state_variables)) {
        if (process %in% names(state_variables[[name]]$rules)) {
          value <- .apply_rule(
            state_variables[[name]],
            name,
            process,
            after[[name]]
          )
          set(after, j = name, value = value)
        }
      }
      after.keys <- as.key(after)
      outside <- which(!after.keys %in% grid.keys)
      if (length(outside) > 0) {
        .stop_outside_states(
          before[outside[1]],
          after[outside[1]],
          process,
          state_variables
        )
      }
      next.keys <- c(next.keys, after.keys)
    }
    new.keys <- setdiff(unique(next.keys), reached)
    reached <- c(reached, new.keys)
    frontier <- grid[grid.keys %in% new.keys]
  }
  depend.matrix[as.key(depend.matrix) %in% reached]
}

# The error for a jump of `process` from state `before` to `after` outside
# the state space, as a state variable's rule sets it to a value it does not
# take in the data up to tau.
.stop_outside_states <- function(before, after, process, state_variables) {
  for (name in names(state_variables)) {
    if (!identical(as.numeric(before[[name]]), as.numeric(after[[name]]))) {
      stop(
        "The rule for ",
        process,
        " of the state variable ",
        name,
        " takes it from ",
        format(before[[name]]),
        " to ",
        format(after[[name]]),
        ", which it does not take in the data up to tau, in a state ",
        "where ",
        process,
        " is at risk. Restrict the rule to the values in the data, ",
        "e.g. with pmin(), or the process with an at_risk condition.",
        call. = FALSE
      )
    }
  }
  stop(
    "Internal error: a jump of ",
    process,
    " leaves the state space.",
    call. = FALSE
  )
}

# The process owning each jump time variable "T.<process>[.<k>]".
.jump_owner <- function(T.vars) {
  vapply(
    strsplit(as.character(T.vars), ".", fixed = TRUE),
    `[`,
    character(1),
    2
  )
}

# The states to report progress for when verbose.
.states_to_report <- function(states) {
  n <- length(states)
  if (n <= 10) {
    return(states)
  }
  floor(seq(1, n, length = n / if (n < 25) 3 else 10))
}

# Predict the hazard increments of the non-censoring processes in every state
# of `depend.matrix`, under baseline treatment `aa` if `intervene.A0`, as
# columns `P[.a<aa>].<process>.<state>` of tmp.long. With HAL fits, also
# sets `state` to the last state each row of tmp.long is in; otherwise
# tmp.long has its `state` already (see .match_states()).
.predict_states <- function(
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
) {
  if (intervene.A0) {
    tmp.long[, (A0.name) := aa]
    tmp.long[,
      (paste0("clever.weight.a", aa)) := clever.weight *
        (A0.obs == aa) /
        (pi.A0.1^(A0.obs == 1) * (1 - pi.A0.1)^(A0.obs == 0))
    ]
    if (!is.null(hal)) {
      hal <- .hal_design_set_arm(hal, tmp.long)
      for (name in intersect(targets, names(which(hal$which.hal)))) {
        tmp.long[,
          (paste0("P.a", aa, ".", name)) := as.numeric(
            .hal_predict(hal, name, tmp.long[["risk.time"]])
          )
        ]
      }
    } else {
      for (name in targets) {
        tmp.long[,
          (paste0("P.a", aa, ".", name)) := .cox_increment(
            cox[[name]][["fit"]],
            tmp.long,
            tmp.long[[paste0("dhazard.", name)]],
            if (name %in% names(at.risks)) {
              tmp.long[[paste0("at.risk.", name)]]
            }
          )
        ]
      }
    }
  }

  tmp.long[, row_id := .I]
  keep <- unique(c(
    "state",
    "id",
    "row_id",
    "risk.time",
    "delta",
    "tstart",
    "time",
    "A0.obs",
    grep("dhazard.", names(tmp.long), value = TRUE, fixed = TRUE),
    grep("at.risk", names(tmp.long), value = TRUE, fixed = TRUE),
    "grid.time",
    varnames,
    hal$N.vars
  ))
  tmp.state <- tmp.long[, names(tmp.long) %in% keep, with = FALSE]

  prefix <- if (intervene.A0) paste0("P.a", aa, ".") else "P."
  dhazards <- lapply(
    stats::setNames(nm = targets),
    function(name) tmp.long[[paste0("dhazard.", name)]]
  )
  state.vars <- setdiff(names(depend.matrix), "state")
  states <- unique(depend.matrix[["state"]])
  # One task per worker, each with a run of consecutive states.
  chunks <- .chunks(states, .n_workers())

  outs <- .parallel_lapply(
    chunks,
    .predict_state_chunk,
    tmp.state = tmp.state,
    observed = if (is.null(hal)) NULL else tmp.long[, state.vars, with = FALSE],
    depend.matrix = depend.matrix,
    state.vars = state.vars,
    targets = targets,
    prefix = prefix,
    cox = cox,
    hal = hal,
    at.risks = at.risks,
    dhazards = dhazards,
    report = if (verbose) .states_to_report(states) else integer(),
    n.states = max(depend.matrix[["state"]])
  )

  for (out in outs) {
    cols <- setdiff(names(out), c("row_id", "state"))
    tmp.long[, (cols) := out[, cols, with = FALSE]]
    if (!is.null(hal)) {
      tmp.long[, state := pmax(state, out[["state"]])]
    }
  }
  tmp.long
}

# The hazard increments of .predict_states() in the states `chunk`, for each
# row of `tmp.state` with its state variables set to those of the state, as
# columns `<prefix><process>.<state>`, and `state`, the last of the states
# each row is in by its `observed` state variables (if given). Progress is
# reported for the states in `report` (only shown when not run on a worker).
.predict_state_chunk <- function(
  chunk,
  tmp.state,
  observed,
  depend.matrix,
  state.vars,
  targets,
  prefix,
  cox,
  hal,
  at.risks,
  dhazards,
  report,
  n.states
) {
  out <- data.table(row_id = tmp.state[["row_id"]], state = 0)
  for (state.jj in chunk) {
    if (state.jj %in% report) {
      message("State ", state.jj, "/", n.states)
    }
    row <- depend.matrix[state == state.jj]
    for (varname in state.vars) {
      tmp.state[[varname]] <- row[[varname]]
    }
    if (!is.null(hal)) {
      hal.state <- .hal_design_set_state(hal, tmp.state)
    }

    for (name in targets) {
      colname <- paste0(prefix, name, ".", state.jj)
      at.risk <- if (name %in% names(at.risks)) at.risks[[name]](row) else TRUE
      if (!is.null(hal) && hal$which.hal[[name]]) {
        P <- as.numeric(
          .hal_predict(hal.state, name, tmp.state[["risk.time"]])
        ) *
          at.risk
      } else {
        P <- .cox_increment(
          cox[[name]][["fit"]],
          tmp.state,
          dhazards[[name]],
          at.risk
        )
      }
      out[[colname]] <- P
    }
    if (!is.null(observed)) {
      in.state <- Reduce(
        `&`,
        lapply(state.vars, function(name) observed[[name]] >= row[[name]]),
        rep(TRUE, nrow(tmp.state))
      )
      out[["state"]][in.state] <- state.jj
    }
  }
  out
}

# Add to the state space one row per value of the time of each jump of
# `specs` (see .jump_specs(); as columns `time.<name>`, Inf when the jump
# has not happened yet), for the states just before that jump.
.add_jump_time_rows <- function(depend.matrix, tmp.long, specs, zero) {
  for (spec in specs) {
    time.col <- spec$time.col
    time.values <- unique(tmp.long[[time.col]])
    if (zero) {
      time.values <- unique(c(0, time.values))
    }
    before.jump <- .jump_count(depend.matrix, spec) == spec$k - 1
    depend.matrix[, (time.col) := Inf]
    depend.matrix[before.jump, (time.col) := time.values[1]]
    if (length(time.values) > 1) {
      base <- depend.matrix[before.jump]
      new.rows <- base[rep(
        seq_len(nrow(base)),
        times = length(time.values) - 1
      )]
      new.rows[, (time.col) := rep(time.values[-1], each = nrow(base))]
      depend.matrix <- rbind(depend.matrix, new.rows)
      depend.matrix <- depend.matrix[order(state, depend.matrix[[time.col]])]
    }
  }
  depend.matrix
}

# Link each state to the state after a jump of each process in
# `process.depend` (columns `s.<process>` of the state space), and each row
# of tmp.long to its row of the state space (`state.row.index`). A jump adds
# one to the count of its process, capped at `count.caps` (by default the
# largest count in the state space). A jump that is the kth event of the
# processes of a jump time variable in `specs` (see .jump_specs()) sets the
# variable to its value at the time of the jump, and a jump sets each of the
# `state_variables` with a rule for its process by the rule. A jump to a
# state not in the state space, which the process is not at risk of, is a
# jump to the same state.
.link_states <- function(
  depend.matrix,
  tmp.long,
  var.depend,
  process.depend,
  specs,
  state_variables = list(),
  count.caps = NULL
) {
  depend.matrix[, state.row.index := seq_len(.N)]
  time.cols <- grep("^time\\.", names(depend.matrix), value = TRUE)
  lookup <- unique(depend.matrix[, c("state", var.depend), with = FALSE])
  for (name in var.depend) {
    set(lookup, j = name, value = as.numeric(lookup[[name]]))
  }
  caps <- lapply(
    stats::setNames(nm = process.depend),
    function(process) {
      if (process %in% names(count.caps)) {
        count.caps[[process]]
      } else {
        max(depend.matrix[[process]])
      }
    }
  )

  for (state.jj in unique(depend.matrix[["state"]])) {
    rows.jj <- depend.matrix[state == state.jj]
    for (process in process.depend) {
      after.jump <- rows.jj[, var.depend, with = FALSE]
      for (name in var.depend) {
        set(after.jump, j = name, value = as.numeric(after.jump[[name]]))
      }
      set(
        after.jump,
        j = process,
        value = pmin(after.jump[[process]] + 1, caps[[process]])
      )
      for (T.var in names(specs)) {
        spec <- specs[[T.var]]
        if (
          process %in%
            spec$processes &&
            all(.jump_count(rows.jj, spec) == spec$k - 1)
        ) {
          set(after.jump, j = T.var, value = rows.jj[[spec$time.col]])
        }
      }
      for (name in names(state_variables)) {
        variable <- state_variables[[name]]
        if (process %in% names(variable$rules)) {
          set(
            after.jump,
            j = name,
            value = .apply_rule(variable, name, process, after.jump[[name]])
          )
        }
      }
      next.state <- lookup[after.jump, on = var.depend, state, mult = "first"]
      next.state[is.na(next.state)] <- state.jj
      depend.matrix[state == state.jj, (paste0("s.", process)) := next.state]
    }

    if (length(time.cols) > 0) {
      timed <- time.cols[vapply(
        time.cols,
        function(time.col) any(rows.jj[[time.col]] < Inf),
        logical(1)
      )]
      if (length(timed) > 0) {
        tmp.long[
          state == state.jj,
          state.row.index := rows.jj[,
            c("state.row.index", timed),
            with = FALSE
          ][tmp.long[state == state.jj], on = timed][["state.row.index"]]
        ]
      } else {
        tmp.long[
          state == state.jj,
          state.row.index := rows.jj[["state.row.index"]]
        ]
      }
    }
  }
  if (length(time.cols) == 0) {
    # One row per state, so a single join links tmp.long to the state space.
    tmp.long[
      depend.matrix,
      on = "state",
      state.row.index := i.state.row.index
    ]
  }

  depend.matrix[, state.row.index := NULL]
  setindex(depend.matrix, NULL)
  list(depend.matrix = depend.matrix, tmp.long = tmp.long)
}
