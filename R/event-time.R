# event_time() terms in the hazard models of fit_initial().
#
# A term such as I(event_time(illness, 1) > 0.5) is a function of the time
# of the first illness. fit_initial() replaces it by the jump time variable
# T.illness, with the term as its function of the time, fitted as a factor.
# Each value of the term is a separate state, so it must take whole numbers
# 0, 1, 2, ... (categories of the time, such as an indicator). Before the
# event, event_time() is 0, so the term is its value at 0.
#
# event_time(c(a, b), k) is the time of the kth event of any of the
# processes a and b, with jump time variable T.a_b (T.a_b.<k> for k > 1).

# The models' right-hand sides with each event_time() term replaced by its
# jump time variable, and the functions of the time defining those
# variables (.fit_initial_impl()'s `depend_time`), with the processes whose
# events each counts (`jump_processes`). Both are named by the variable
# without its "T." prefix.
# `times` are the event times in the data, to check the terms at, and
# `columns` the data columns. Returns list(formulas, depend_time,
# jump_processes), `formulas` named by process.
.event_time_covariates <- function(models, times, columns, max_count) {
  found <- new.env()
  found$depend_time <- list()
  found$jump_processes <- list()
  found$exprs <- list()
  found$terms <- list()
  check_times <- sort(unique(c(0, times)))
  formulas <- lapply(models, function(model) {
    formula <- model$formula
    formula[[3]] <- .replace_event_time(
      formula[[3]],
      found,
      env = environment(formula),
      process_names = names(models),
      columns = columns,
      check_times = check_times,
      max_count = max_count,
      process = model$name
    )
    formula
  })
  list(
    formulas = formulas,
    depend_time = found$depend_time,
    jump_processes = found$jump_processes
  )
}

# `expr`, a formula's right-hand side, with each term containing
# event_time() replaced by its jump time variable, which is recorded in
# `found`. As in .replace_terms(), looks into the formula operators that
# combine terms but not into other calls.
.replace_event_time <- function(
  expr,
  found,
  env,
  process_names,
  columns,
  check_times,
  max_count,
  process
) {
  operators <- c("+", "-", ":", "*", "/", "^", "%in%", "(")
  if (!is.call(expr)) {
    return(expr)
  }
  if (as.character(expr[[1]])[1] %in% operators) {
    for (i in seq_along(expr)[-1]) {
      expr[[i]] <- .replace_event_time(
        expr[[i]],
        found,
        env,
        process_names,
        columns,
        check_times,
        max_count,
        process
      )
    }
    return(expr)
  }
  calls <- .event_time_calls(expr)
  if (length(calls) == 0) {
    return(expr)
  }

  term <- .deparse_one_line(expr)
  what <- paste0("The term ", term, " in the model of ", process)
  if (length(unique(vapply(calls, .deparse_one_line, character(1)))) > 1) {
    stop(what, " uses several event_time()s; use one per term.", call. = FALSE)
  }
  jump <- .parse_event_time(calls[[1]], process_names, env, what)
  if (jump$k > max_count) {
    stop(
      what,
      " uses event ",
      jump$k,
      " of ",
      jump$label,
      ", but only max_count = ",
      max_count,
      " events are tracked; increase max_count.",
      call. = FALSE
    )
  }
  indicator_hint <- paste0(
    "; it must take whole numbers 0, 1, 2, ... (categories of the time), ",
    "such as I(event_time(",
    jump$arg,
    ", ",
    jump$k,
    ") > 0.5), as each value is a separate state in the estimation."
  )
  if (identical(expr, calls[[1]])) {
    stop(what, " is the event time itself", indicator_hint, call. = FALSE)
  }

  # The term as a function of the event time, `..time..`.
  g_expr <- .substitute_call(expr, calls[[1]], as.name("..time.."))
  others <- setdiff(all.vars(g_expr), "..time..")
  bad <- others[
    others %in%
      c(columns, process_names, "time", "t") |
      !vapply(others, exists, logical(1), envir = env)
  ]
  if (length(bad) > 0) {
    stop(
      what,
      " uses ",
      paste(bad, collapse = ", "),
      ", but besides event_time() it may only use constants.",
      call. = FALSE
    )
  }
  g <- function(time) {
    values <- eval(g_expr, list(..time.. = time), env)
    if (length(values) == 1) {
      values <- rep(values, length(time))
    }
    as.numeric(values)
  }
  values <- tryCatch(g(check_times), error = function(e) {
    stop(what, " could not be evaluated: ", conditionMessage(e), call. = FALSE)
  })
  if (length(values) != length(check_times) || anyNA(values)) {
    stop(what, " is not one value per time", indicator_hint, call. = FALSE)
  }
  if (!checkmate::test_integerish(values, lower = 0)) {
    stop(
      what,
      " takes values other than whole numbers 0, 1, 2, ...",
      indicator_hint,
      call. = FALSE
    )
  }

  base <- paste(jump$processes, collapse = "_")
  key <- if (jump$k == 1) base else paste0(base, ".", jump$k)
  g_string <- .deparse_one_line(g_expr)
  if (key %in% names(found$exprs) && found$exprs[[key]] != g_string) {
    stop(
      "The models use several functions of event ",
      jump$k,
      " of ",
      jump$label,
      " (",
      found$terms[[key]],
      " and ",
      term,
      "); only one per event is supported.",
      call. = FALSE
    )
  }
  found$exprs[[key]] <- g_string
  found$terms[[key]] <- term
  found$depend_time[[key]] <- g
  found$jump_processes[[key]] <- jump$processes
  as.name(paste0("T.", key))
}

# The event_time() calls in `expr`.
.event_time_calls <- function(expr) {
  if (!is.call(expr)) {
    return(list())
  }
  if (identical(expr[[1]], as.name("event_time"))) {
    return(list(expr))
  }
  unlist(
    lapply(as.list(expr)[-1], .event_time_calls),
    recursive = FALSE
  )
}

# The processes and event number of an event_time(proc, k) call, where
# `proc` is a process or c() of processes, and `k` defaults to 1. Returns
# list(processes, k, arg, label): the processes in the order of
# `process_names`, `proc` as written, and a description of the processes.
.parse_event_time <- function(call, process_names, env, what) {
  call <- match.call(function(proc, k = 1) NULL, call)
  proc <- call$proc
  arg <- .deparse_one_line(proc)
  names <- if (is.name(proc)) {
    as.character(proc)
  } else if (
    is.call(proc) &&
      identical(proc[[1]], as.name("c")) &&
      length(proc) > 1 &&
      all(vapply(as.list(proc)[-1], is.name, logical(1)))
  ) {
    vapply(as.list(proc)[-1], as.character, character(1))
  } else {
    NULL
  }
  if (is.null(names)) {
    stop(
      what,
      " uses event_time() of ",
      arg,
      "; give a process, or several as c(<process>, <process>).",
      call. = FALSE
    )
  }
  unknown <- setdiff(names, process_names)
  if (length(unknown) > 0) {
    stop(
      what,
      " uses event_time() of ",
      paste(unknown, collapse = ", "),
      ", which is not a process in `models`.",
      call. = FALSE
    )
  }
  processes <- intersect(process_names, names)
  k <- if (is.null(call$k)) 1 else eval(call$k, env)
  if (!checkmate::test_count(k, positive = TRUE)) {
    stop(
      what,
      " uses event_time() with k = ",
      .deparse_one_line(call$k),
      ", but k must be a positive whole number (the kth event).",
      call. = FALSE
    )
  }
  label <- if (length(processes) == 1) {
    processes
  } else {
    paste("any of", paste(processes, collapse = ", "))
  }
  list(processes = processes, k = as.integer(k), arg = arg, label = label)
}

# `expr` with each subexpression identical to `target` replaced by
# `replacement`.
.substitute_call <- function(expr, target, replacement) {
  if (identical(expr, target)) {
    return(replacement)
  }
  if (is.call(expr)) {
    for (i in seq_along(expr)) {
      expr[[i]] <- .substitute_call(expr[[i]], target, replacement)
    }
  }
  expr
}
