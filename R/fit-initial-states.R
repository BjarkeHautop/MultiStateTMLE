# The discrete state space of fit_initial(): the product of the values of the
# state variables, the hazards predicted in each state, and the links from
# each state to the state after a jump.

# All combinations of the values of the state variables, with process counts
# truncated at max_count. A jump time variable T.<name> takes the values it
# has in the data and, from `time.<name>`, the values it would get from a
# jump at any time of the grid.
.state_grid <- function(tmp.long, state.names, process.names, max_count) {
  grid.list <- lapply(state.names, function(varname) {
    time.col <- sub("^T\\.", "time.", varname)
    if (startsWith(varname, "T.") && time.col %in% names(tmp.long)) {
      return(sort(unique(c(tmp.long[[varname]], tmp.long[[time.col]]))))
    }
    0:min(
      if (varname %in% process.names) max_count else Inf,
      max(tmp.long[[varname]])
    )
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

# Drop the states in which no process with an at-risk condition is at risk,
# or would be after a jump, ignoring conditions that hold in every state.
.prune_states <- function(depend.matrix, at.risks) {
  at.risk <- lapply(at.risks, function(at_risk) {
    at_risk(depend.matrix) | at_risk(depend.matrix - 1)
  })
  at.risk <- Filter(function(x) !all(x), at.risk)
  if (length(at.risk) > 0) {
    depend.matrix <- depend.matrix[Reduce(`|`, at.risk)]
  }
  depend.matrix <- depend.matrix[order(state)]
  depend.matrix[, state := seq_len(.N)]
  depend.matrix
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
# columns `P[.a<aa>].<process>.<state>` of tmp.long. Also sets `state` to
# the last state each row of tmp.long is in.
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
            tmp.long[[paste0("dhazard.", name)]]
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
    observed = tmp.long[, state.vars, with = FALSE],
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
    tmp.long[, state := pmax(state, out[["state"]])]
  }
  tmp.long
}

# The hazard increments of .predict_states() in the states `chunk`, for each
# row of `tmp.state` with its state variables set to those of the state, as
# columns `<prefix><process>.<state>`, and `state`, the last of the states
# each row is in by its `observed` state variables. Progress is reported for
# the states in `report` (only shown when not run on a worker).
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
    in.state <- rep(TRUE, nrow(tmp.state))
    for (varname in state.vars) {
      tmp.state[[varname]] <- row[[varname]]
      in.state <- in.state & (observed[[varname]] >= row[[varname]])
    }
    if (!is.null(hal)) {
      hal.state <- .hal_design_set_state(hal, tmp.state)
    }

    for (name in targets) {
      colname <- paste0(prefix, name, ".", state.jj)
      if (!is.null(hal) && hal$which.hal[[name]]) {
        P <- as.numeric(
          .hal_predict(hal.state, name, tmp.state[["risk.time"]])
        )
      } else {
        P <- .cox_increment(
          cox[[name]][["fit"]],
          tmp.state,
          dhazards[[name]]
        )
      }
      if (name %in% names(at.risks)) {
        P <- P * at.risks[[name]](row)
      }
      out[[colname]] <- P
    }
    out[["state"]][in.state] <- state.jj
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
# of tmp.long to its row of the state space (`state.row.index`). A jump that
# is the kth event of the processes of a jump time variable in `specs` (see
# .jump_specs()) sets the variable to its value at the time of the jump.
.link_states <- function(
  depend.matrix,
  tmp.long,
  var.depend,
  process.depend,
  specs
) {
  depend.matrix[, state.row.index := seq_len(.N)]
  time.cols <- grep("^time\\.", names(depend.matrix), value = TRUE)
  lookup <- unique(depend.matrix[, c("state", var.depend), with = FALSE])
  max.values <- lapply(
    stats::setNames(nm = process.depend),
    function(process) max(depend.matrix[[process]])
  )

  for (state.jj in unique(depend.matrix[["state"]])) {
    rows.jj <- depend.matrix[state == state.jj]
    for (process in process.depend) {
      count <- rows.jj[[process]]
      after.jump <- rows.jj[, c("state.row.index", var.depend), with = FALSE]
      after.jump[[process]] <- pmin(count + 1L, max.values[[process]])
      for (T.var in names(specs)) {
        spec <- specs[[T.var]]
        if (
          process %in%
            spec$processes &&
            all(.jump_count(rows.jj, spec) == spec$k - 1)
        ) {
          after.jump[[T.var]] <- rows.jj[[spec$time.col]]
        }
      }
      next.state <- unique(merge(after.jump, lookup, by = var.depend))[
        order(state.row.index)
      ][["state"]]
      if (length(next.state) == 0) {
        next.state <- state.jj
      }
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
