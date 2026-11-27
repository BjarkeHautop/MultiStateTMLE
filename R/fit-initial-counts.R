# The number of events of each process that the state space of
# fit_initial() tracks.
#
# A process in a model formula or a condition is its number of events so
# far, used as written. The state space is finite, so it tracks each count
# up to a cap: the smallest count from which on nothing that uses the
# process changes. The cap is found by evaluating each use on the data with
# the count set to larger and larger values. A use that keeps changing,
# such as a plain count or a difference of counts, can't be tracked and is
# an error.

# How far beyond the largest count in the data, or the largest count a use
# needs, a use must be seen not to change.
.count_margin <- 10

# The range of counts to check each process at: its limit, 1 for a process
# that jumps at most once in the data, and 0 for a terminal one. Inf for a
# recurrent process without a limit. A named numeric vector.
.count_upper <- function(processes, process.types, limits) {
  upper <- vapply(
    processes,
    function(process) {
      limit <- limits[[process]]
      if (!is.null(limit)) {
        return(as.numeric(limit))
      }
      switch(
        process.types[[process]],
        terminal = 0,
        one_jump = 1,
        Inf
      )
    },
    numeric(1)
  )
  stats::setNames(upper, processes)
}

# The uses of the processes in the hazard models and their at-risk
# conditions, as a list of list(what, fun, processes): a description, a
# function of the data, and the processes it reads. `models` are the
# one-sided formulas of the covariates, named by process.
.count_uses <- function(models, at.risks, processes) {
  uses <- list()
  for (name in names(models)) {
    model <- models[[name]]
    variables <- as.list(attr(stats::terms(model), "variables"))[-1]
    for (variable in variables) {
      reads <- intersect(all.vars(variable), processes)
      if (length(reads) == 0) {
        next
      }
      uses[[length(uses) + 1]] <- list(
        what = paste0(
          "The term ",
          .deparse_one_line(variable),
          " in the model of ",
          name
        ),
        fun = .eval_function(variable, environment(model)),
        processes = reads
      )
    }
  }
  for (name in names(at.risks)) {
    at_risk <- at.risks[[name]]
    reads <- intersect(.at_risk_columns(at_risk), processes)
    if (length(reads) > 0) {
      uses[[length(uses) + 1]] <- list(
        what = paste0("The at-risk condition of ", name),
        fun = at_risk,
        processes = reads
      )
    }
  }
  uses
}

# A function of a data list or data.table evaluating `expr` in it, with
# other variables from `env`.
.eval_function <- function(expr, env) {
  force(expr)
  force(env)
  function(data) eval(expr, as.list(data), env)
}

# The values of a use, comparable with identical(): factors as their
# labels, and without attributes.
.comparable <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  attributes(x) <- NULL
  x
}

# The smallest count `c` of `process` such that `fun` gives the same
# values on `data` with the count set to any of c, ..., `upper`.
.constant_from <- function(fun, data, process, upper) {
  data <- as.list(data)
  n <- length(data[[process]])
  values <- lapply(0:upper, function(count) {
    data[[process]] <- rep(count, n)
    .comparable(suppressWarnings(fun(data)))
  })
  from <- upper
  while (from > 0 && identical(values[[from]], values[[upper + 1]])) {
    from <- from - 1
  }
  from
}

# The cap of each of `processes`: the smallest count from which on none of
# `uses` (see .count_uses()) changes on `dt`, and at least `needed`. `upper`
# (see .count_upper()) bounds the counts; for a recurrent process without a
# limit, the uses are checked up to .count_margin beyond the largest count
# in the data (`observed`). A use of several processes is checked for each
# count of the others in their whole range, as the state space holds every
# combination. A named integer vector.
.count_caps <- function(dt, processes, uses, upper, observed, needed) {
  open <- is.infinite(upper[processes])
  tops <- ifelse(
    open,
    pmax(observed[processes], needed[processes]) + .count_margin,
    upper[processes]
  )
  names(open) <- names(tops) <- processes
  data <- as.list(dt)
  n <- length(data[["id"]])

  caps <- stats::setNames(integer(length(processes)), processes)
  for (process in processes) {
    top <- tops[[process]]
    cap <- min(needed[[process]], top)
    for (use in uses) {
      if (!process %in% use$processes) {
        next
      }
      others <- setdiff(use$processes, process)
      combos <- expand.grid(lapply(
        stats::setNames(nm = others),
        function(other) 0:tops[[other]]
      ))
      for (i in seq_len(max(1, nrow(combos)))) {
        for (other in others) {
          data[[other]] <- rep(combos[[other]][i], n)
        }
        from <- .constant_from(use$fun, data, process, top)
        if (open[[process]] && from > top - .count_margin) {
          stop(.keeps_changing(use, process), call. = FALSE)
        }
        cap <- max(cap, from)
      }
      for (other in others) {
        data[[other]] <- dt[[other]]
      }
    }
    caps[[process]] <- as.integer(cap)
  }
  caps
}

# The error message for a use that keeps changing with the number of events
# of `process`.
.keeps_changing <- function(use, process) {
  paste0(
    use$what,
    " keeps changing with the number of ",
    process,
    " events, but the state space can only track a bounded number of them. ",
    if (length(use$processes) > 1) {
      paste0(
        "For a quantity of several processes, such as a difference of ",
        "counts, define a state_variable(). "
      )
    },
    "Use a term that stops changing, such as I(",
    process,
    " >= 1), pmin(",
    process,
    ", 3) or factor(pmin(",
    process,
    ", 3)) for an effect of each count up to 3."
  )
}

# Check that `fun` of `data` (with the processes' counts in `caps` capped)
# does not change when a count goes beyond its cap, up to `upper` (see
# .count_upper()) or .count_margin beyond the cap. `what` describes `fun`
# in the error message.
.check_beyond_caps <- function(fun, data, caps, upper, what) {
  data <- as.list(data)
  for (process in names(caps)) {
    if (is.null(data[[process]])) {
      next
    }
    cap <- caps[[process]]
    top <- min(upper[[process]], cap + .count_margin)
    if (top <= cap) {
      next
    }
    n <- length(data[[process]])
    at <- function(count) {
      data[[process]] <- rep(count, n)
      .comparable(suppressWarnings(fun(data)))
    }
    at.cap <- at(cap)
    for (count in (cap + 1):top) {
      if (!identical(at(count), at.cap)) {
        stop(
          what,
          " changes with the number of ",
          process,
          " events beyond ",
          cap,
          ", but the state space only tracks ",
          cap,
          " (the most that the hazard models tell apart). Use a ",
          "condition that does not change beyond ",
          cap,
          " events, or a state_variable().",
          call. = FALSE
        )
      }
    }
  }
  invisible(TRUE)
}

# The largest number of events of each of `processes` (with delta codes
# `deltas`) that an individual has in the data `dt`. A named numeric vector.
.max_counts <- function(dt, processes, deltas) {
  counts <- vapply(
    deltas,
    function(d) {
      max(c(0, dt[, list(n = sum(delta == d)), by = "id"][["n"]]))
    },
    numeric(1)
  )
  stats::setNames(counts, processes)
}

# The number of events of each process that the jump time variables of
# `specs` (see .jump_specs()) need tracked: for the kth event of several
# processes, k of each, or as many as it can have by `upper` (see
# .count_upper()). It is an error if the processes can't have k events
# together. A named numeric vector for the processes of `upper`.
.needed_counts <- function(specs, upper) {
  needed <- stats::setNames(numeric(length(upper)), names(upper))
  for (spec in specs) {
    most <- sum(upper[spec$processes])
    if (most < spec$k) {
      stop(
        "An event_time() term uses event ",
        spec$k,
        " of ",
        paste(spec$processes, collapse = " and "),
        ", which ",
        if (length(spec$processes) > 1) "have" else "has",
        " at most ",
        most,
        " (by the limits, or in the data).",
        call. = FALSE
      )
    }
    for (process in spec$processes) {
      needed[[process]] <- max(
        needed[[process]],
        min(spec$k, upper[[process]])
      )
    }
  }
  needed
}
