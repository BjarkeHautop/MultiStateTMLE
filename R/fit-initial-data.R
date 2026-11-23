# Data preparation for fit_initial(): the variables and types of the
# processes, the history columns, and the expansion of the data to the grid
# of event times.

# The variables used by the hazard models and the at-risk conditions: the
# variables, not the terms, so that e.g. log(biomarker) gives biomarker.
.model_variables <- function(fit_types, at.risks) {
  model.vars <- lapply(fit_types, function(fit.type) {
    rhs <- strsplit(fit.type[["model"]], "~")[[1]][2]
    all.vars(as.formula(paste("~", rhs)))
  })
  unique(unlist(
    c(model.vars, lapply(at.risks, .at_risk_columns)),
    use.names = FALSE
  ))
}

# The columns an at-risk function reads, as `dt[["<name>"]]`.
.at_risk_columns <- function(at_risk_fun) {
  cols <- character()
  recurse <- function(e) {
    if (
      is.call(e) &&
        identical(e[[1]], as.name("[[")) &&
        length(e) == 3 &&
        is.character(e[[3]])
    ) {
      cols <<- c(cols, e[[3]])
    }
    if (is.call(e) || is.expression(e)) {
      lapply(as.list(e), recurse)
    }
  }
  recurse(body(at_risk_fun))
  unique(cols)
}

# The type of each process, inferred from the data: "terminal" if it only
# occurs as the last event of an individual, "recurrent" if it occurs more
# than once for some individual, and "one_jump" otherwise. A named list.
.process_types <- function(dt, process.deltas, process.names) {
  by.id <- dt[, list(delta, before.end = time < max(time)), by = "id"]
  types <- vapply(
    process.deltas,
    function(d) {
      if (!any(by.id[["delta"]] == d & by.id[["before.end"]])) {
        "terminal"
      } else if (any(dt[delta == d, .N, by = "id"][["N"]] > 1)) {
        "recurrent"
      } else {
        "one_jump"
      }
    },
    character(1)
  )
  stats::setNames(as.list(types), process.names)
}

# The jump number of a state or jump time variable such as "z.2" (or
# "T.z.2"), read from dot-separated part `position`; 1 if it has none.
.jump_number <- function(name, position) {
  value <- suppressWarnings(
    as.numeric(strsplit(name, ".", fixed = TRUE)[[1]][position])
  )
  if (is.na(value)) 1 else value
}

# Add to `x` (sorted by id and time, with `delta`) the number of jumps of
# each process in `processes` before each row.
.add_count_columns <- function(x, processes, deltas) {
  for (j in seq_along(processes)) {
    process <- processes[j]
    d <- deltas[j]
    x[, (process) := c(0L, cumsum(delta == d)[-.N]), by = "id"]
  }
  x
}

# Add to `x` the jump time column `T.<name>` of each state variable in
# `depend_time` (`<process>` or `<process>.<k>`, the kth jump): the
# depend_time function at 0 before the jump, as event_time() is 0 then, and
# at the time of the jump after it. Whether the jump has happened is known
# from the count of the process, so the value before it can be any.
# `deltas` are the delta codes of the processes, named by process. With
# `grid = TRUE`, also add `time.<name>`, the function at the time of each
# row.
.add_jump_time_columns <- function(x, depend_time, deltas, grid = FALSE) {
  for (time.var in names(depend_time)) {
    time.fun <- depend_time[[time.var]]
    T.var <- paste0("T.", time.var)
    d <- deltas[[strsplit(time.var, ".", fixed = TRUE)[[1]][1]]]
    k <- .jump_number(time.var, 2)
    x[,
      (T.var) := {
        jumps <- time[delta == d]
        if (length(jumps) >= k) jumps[k] else Inf
      },
      by = "id"
    ]
    before <- x[["tstart"]] < x[[T.var]]
    x[(before), (T.var) := as.numeric(time.fun(0))]
    x[(!before), (T.var) := as.numeric(time.fun(get(T.var)))]
    if (grid) {
      x[, (paste0("time.", time.var)) := as.numeric(time.fun(time))]
    }
  }
  x
}

# Expand the data `dt` to one row per individual and time in `unique.times`,
# carrying the covariates `varnames` forward, and join the baseline hazard
# increments `dhazard.<process>` of the Cox fits.
.expand_to_grid <- function(dt, hazards, varnames, unique.times) {
  all.times <- data.table(expand.grid(
    time = unique.times,
    id = unique(dt[["id"]])
  ))
  tmp.long <- merge(
    dt[, time.obs := time],
    all.times,
    by = c("id", "time"),
    all = TRUE
  )
  setorder(tmp.long, id, time.obs, na.last = TRUE)

  for (varname in varnames) {
    tmp.long[, (varname) := nafill(get(varname), type = "locf"), by = "id"]
    if (anyNA(tmp.long[[varname]])) {
      stop(
        "Unexpected leading NA in '",
        varname,
        "' after forward-fill (first observation missing for some id)."
      )
    }
  }

  setorder(tmp.long, id, time, na.last = TRUE)
  tmp.long[is.na(delta), delta := 0]
  tmp.long[, tstop := NULL]

  # Join the (small) per-process hazard tables together first, and then onto
  # tmp.long in a single pass.
  combined.hazard <- Reduce(
    function(a, b) merge(a, b, by = "time", all = TRUE),
    hazards
  )
  tmp.long <- merge(tmp.long, combined.hazard, by = "time", all = TRUE)
  setorder(tmp.long, id, time, na.last = TRUE)
  tmp.long <- tmp.long[!is.na(id)]

  tmp.long[, time.obs := nafill(time.obs, "nocb"), by = "id"]
  tmp.long[is.na(time.obs), time.obs := -Inf]
  tmp.long[, final.time := max(time.obs), by = "id"]

  for (col in grep("^dhazard\\.", names(tmp.long), value = TRUE)) {
    set(tmp.long, which(is.na(tmp.long[[col]])), col, 0)
  }
  tmp.long
}
