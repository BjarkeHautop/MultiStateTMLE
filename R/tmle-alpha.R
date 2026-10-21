#' Targeted estimation under a stochastic intensity intervention
#'
#' Estimates a target parameter under a stochastic intensity intervention
#' that scales the intensity of the process `z_name` by `alpha`, optionally
#' combined with an intervention `a` on the baseline treatment. The initial
#' fit from [fit_initial()] is updated by iterative targeted maximum
#' likelihood estimation (TMLE).
#'
#' By default the target parameter is the risk of a `target` event by the
#' horizon `tau` given to [fit_initial()] (or the expected number of events
#' for a recurrent `target`).
#'
#' @param initial_fit An initial fit, as returned by [fit_initial()].
#' @param target String. Name of the target process.
#' @param alpha Non-negative number. Multiplicative intervention on the
#'   intensity of `z_name`. `alpha = 1` is no intervention. Ignored if
#'   `alpha_list` is given.
#' @param z_name String. Name of the process whose intensity `alpha` scales.
#'   If no such process exists in `initial_fit`, `alpha` is ignored with a
#'   warning.
#' @param alpha_list Experimental. Optional named list of functions, one per
#'   process, giving a time-varying intervention on that process's intensity.
#'   Each function's arguments name the columns of `initial_fit$tmp.long` it
#'   depends on. Supersedes `alpha` and `z_name`.
#' @param a Optional baseline treatment value (`0` or `1`) to intervene on.
#'   Must be one of the values of `a` given to [fit_initial()].
#' @param years_lost Optional positive integer. If given, the target
#'   parameter is the expected time lost (the integral of the risk of
#'   `target` up to `tau`) instead of the risk at `tau`. The integral is
#'   approximated with blocks of `years_lost` time points.
#' @param only_first Logical. If `TRUE`, only the first `target` event
#'   counts, and only the first jump of `z_name` is intervened on.
#' @param target_only_in_state Optional function of the state space
#'   (`initial_fit$depend.matrix`) returning a logical vector, marking the
#'   states in which `target` events count.
#' @param target_by_state Logical. If `TRUE`, use state-specific clever
#'   covariates in the targeting step.
#' @param conv_const Positive number. The TMLE iterations stop when the
#'   absolute mean of the efficient influence curve falls below
#'   `conv_const * se / log(n)`.
#' @param one_step Logical. If `TRUE`, return the one-step estimator instead
#'   of iterating the TMLE update.
#' @param max_iter Positive integer. Maximum number of TMLE iterations.
#' @param min_iter Non-negative integer. Minimum number of TMLE iterations
#'   before checking convergence.
#' @param truncate_weights Non-negative number. If positive, the clever
#'   weights of the intervention are truncated at this value.
#' @param output_convergence Logical. If `TRUE`, include convergence
#'   diagnostics in the output.
#' @param output_eic Logical. If `TRUE`, include the efficient influence
#'   curve in the output.
#' @param output_weights Optional numeric vector of probabilities at which to
#'   report quantiles of the per-individual maximal intervention and
#'   censoring weights.
#' @param output_a_weights Optional numeric vector of probabilities at which
#'   to report quantiles of the per-individual maximal treatment weights.
#' @param cores Positive integer. Number of cores used by
#'   [parallel::mclapply()] for the per-individual computations.
#' @param verbose Logical. If `TRUE`, print progress.
#'
#' @return A list with element `estimate`, a named vector with the estimate
#'   (`tmle.est`, or `one.step.est` if `one_step = TRUE`), its standard
#'   error `se`, and the initial plug-in estimate `g.est`. Depending on the
#'   arguments, the list also contains `convergence`, `weights`,
#'   `cens.weights`, `a.weights` and `eic`.
#'
#' @examples
#' fit <- fit_initial(
#'   illness_death,
#'   fit_types = list(
#'     z = list(
#'       model = "Surv(tstart, tstop, delta == 2) ~ L0",
#'       at_risk = function(dt) dt[["z"]] == 0
#'     ),
#'     death = "Surv(tstart, tstop, delta == 1) ~ L0 + z",
#'     censoring = "Surv(tstart, tstop, delta == 0) ~ L0"
#'   ),
#'   tau = 1
#' )
#'
#' # Risk of death by tau = 1 if the intensity of z were halved.
#' est <- tmle_alpha(fit, target = "death", alpha = 0.5, z_name = "z")
#' est$estimate
#' @export
tmle_alpha <- function(
  initial_fit,
  target,
  alpha = 1,
  z_name = "z",
  alpha_list = NULL,
  a = NULL,
  years_lost = NULL,
  only_first = FALSE,
  target_only_in_state = NULL,
  target_by_state = FALSE,
  conv_const = 1,
  one_step = FALSE,
  max_iter = 10,
  min_iter = 0,
  truncate_weights = 0,
  output_convergence = FALSE,
  output_eic = FALSE,
  output_weights = NULL,
  output_a_weights = NULL,
  cores = 1,
  verbose = FALSE
) {
  .assert_initial_fit(initial_fit)
  checkmate::assert_string(target)
  checkmate::assert_number(alpha, lower = 0, finite = TRUE)
  checkmate::assert_string(z_name)
  checkmate::assert_list(
    alpha_list,
    types = "function",
    names = "unique",
    null.ok = TRUE
  )
  checkmate::assert_int(a, lower = 0, upper = 1, null.ok = TRUE)
  checkmate::assert_count(years_lost, positive = TRUE, null.ok = TRUE)
  checkmate::assert_flag(only_first)
  checkmate::assert_function(target_only_in_state, null.ok = TRUE)
  checkmate::assert_flag(target_by_state)
  checkmate::assert_number(conv_const, lower = 0, finite = TRUE)
  checkmate::assert_flag(one_step)
  checkmate::assert_count(max_iter, positive = TRUE)
  checkmate::assert_count(min_iter)
  checkmate::assert_number(truncate_weights, lower = 0, finite = TRUE)
  checkmate::assert_flag(output_convergence)
  checkmate::assert_flag(output_eic)
  checkmate::assert_numeric(
    output_weights,
    lower = 0,
    upper = 1,
    any.missing = FALSE,
    null.ok = TRUE
  )
  checkmate::assert_numeric(
    output_a_weights,
    lower = 0,
    upper = 1,
    any.missing = FALSE,
    null.ok = TRUE
  )
  checkmate::assert_count(cores, positive = TRUE)
  checkmate::assert_flag(verbose)

  tmp.long <- copy(initial_fit$tmp.long)
  depend.matrix <- copy(initial_fit$depend.matrix)

  n <- length(unique(tmp.long[["id"]]))

  process.names <- initial_fit$process.names
  process.deltas <- initial_fit$process.deltas
  process.types <- initial_fit$process.types

  at.risks <- initial_fit$at.risks

  target.id <- which(tolower(process.names) == tolower(target))
  if (length(target.id) != 1) {
    stop(
      "target must be one of the processes in initial_fit: ",
      paste(process.names, collapse = ", "),
      "."
    )
  }
  target.name <- process.names[target.id]

  z.name <- z_name
  z.process.id <- which(tolower(process.names) == tolower(z_name))

  if (length(alpha_list) > 0) {
    z.present <- TRUE
  } else if (length(z.process.id) == 0) {
    z.present <- FALSE
    z.process.id <- Inf
    z.name <- "not.applied"
    if (!identical(alpha, 1)) {
      warning(
        "No process '",
        z_name,
        "' found in initial_fit; alpha is ignored (treated as 1)."
      )
    }
  } else {
    z.present <- TRUE
    z.name <- process.names[z.process.id]
  }

  #--------------------------------
  #-- initialization step:

  states <- copy(depend.matrix)
  state.names <- names(states)
  setkeyv(states, c("state", state.names[state.names %in% process.names]))

  if (only_first) {
    process.types[[target.name]] <- "one.jump"
    tmp.long <- tmp.long[tmp.long[[target.name]] <= 1]
  }

  tmp.long[, at.risk := (time <= final.time)]

  if (target.name %in% names(at.risks)) {
    if (process.types[target.name] == "one.jump") {
      tmp.long[, at.risk := at.risk * (at.risks[[target.name]](tmp.long))]
    }
    states[, at.risk := at.risks[[target.name]](depend.matrix)]
  }

  if (only_first && z.present) {
    tmp.long[, at.risk := at.risk * (tmp.long[[z.name]] <= 1)]
    states[, at.risk := at.risk & (states[[z.name]] <= 1)]
  }

  if (!is.null(target_only_in_state)) {
    states[, target.only.in.state := target_only_in_state(states)]
  }

  for (process.name in process.names[!(process.names %in% names(at.risks))]) {
    tmp.long[, (paste0("at.risk.", process.name)) := at.risk]
  }

  P.prefix <- if (length(a) > 0) paste0("P.a", a, ".") else "P."

  if (length(alpha_list) > 0) {
    for (alpha.kk in seq_along(alpha_list)) {
      tmp.long[,
        paste0("alpha.", names(alpha_list)[alpha.kk]) := do.call(
          alpha_list[[alpha.kk]],
          .SD
        ),
        .SDcols = names(formals(alpha_list[[alpha.kk]]))
      ]
    }
  }

  #--------------------------------
  #-- compute clever weights:

  if (z.present) {
    if (length(alpha_list) > 0) {
      .update_alpha_list_weights(tmp.long, alpha_list)
    } else {
      .update_alpha_weights(tmp.long, alpha, z.name)
    }
    if (truncate_weights > 0) {
      no.truncated <- tmp.long[
        time <= final.time & clever.weight.alpha > truncate_weights,
        length(unique(id))
      ]
      tmp.long[
        time <= final.time & clever.weight.alpha > truncate_weights,
        clever.weight.alpha := truncate_weights
      ]
    }
  } else {
    tmp.long[, clever.weight.alpha := 1]
  }

  if (length(a) > 0) {
    a.weight.col <- paste0("clever.weight.a", a)
    if (!(a.weight.col %in% names(tmp.long))) {
      stop(
        "Column '",
        a.weight.col,
        "' not found in initial_fit$tmp.long; a must be one of the values ",
        "of a given to fit_initial()."
      )
    }
    tmp.long[, clever.weight := tmp.long[[a.weight.col]]]
  }

  quantiles <- function(x, probs) {
    out <- vapply(probs, function(p) as.numeric(quantile(x, p)), numeric(1))
    names(out) <- paste0("q", probs * 100)
    out
  }

  if (length(output_weights) > 0) {
    out.weights <- quantiles(
      tmp.long[time <= final.time, max(clever.weight.alpha), by = "id"][[2]],
      output_weights
    )
    if (truncate_weights > 0) {
      out.weights <- c(out.weights, no.truncated = no.truncated)
    }
    out.cens.weights <- quantiles(
      tmp.long[time <= final.time, max(clever.weight), by = "id"][[2]],
      output_weights
    )
  }

  if (length(output_a_weights) > 0) {
    out.a.weights <- quantiles(
      tmp.long[time <= final.time, max(clever.weight), by = "id"][[2]],
      output_a_weights
    )
  }

  #--------------------------------
  #-- counterfactual hazards of the intervened processes:

  intervened <- if (length(alpha_list) > 0) names(alpha_list) else z.name
  for (name.kk in intervened) {
    for (varname in grep(
      paste0(P.prefix, name.kk),
      names(tmp.long),
      value = TRUE
    )) {
      if (varname != paste0(P.prefix, name.kk)) {
        if (length(alpha_list) > 0) {
          tmp.long[, (varname) := get(paste0("alpha.", name.kk)) * get(varname)]
        } else {
          tmp.long[, (varname) := alpha * get(varname)]
        }
      }
    }
  }

  #--------------------------------
  #-- tmle iterations:

  converged <- FALSE
  clever.Q.label <- if (!is.null(years_lost)) {
    "clever.Q.years.lost."
  } else {
    "clever.Q."
  }

  # The intervention's factor on the score of process number jj: alpha for
  # z itself, the alpha_list function for the processes it covers, else 1.
  alpha.factor <- function(jj) {
    name.jj <- process.names[jj]
    if (length(alpha_list) > 0) {
      if (name.jj %in% names(alpha_list)) {
        return(tmp.long[[paste0("alpha.", name.jj)]])
      }
      return(1)
    }
    alpha^(jj == z.process.id)
  }

  # Per-individual score of process jj, with the hazards perturbed by eps
  # (along the clever covariate if target_by_state).
  score <- function(jj, eps = 0) {
    name.jj <- process.names[jj]
    clever <- tmp.long[[paste0(clever.Q.label, name.jj)]]
    fluct <- if (target_by_state) exp(eps * clever) else exp(eps)
    s <- alpha.factor(jj) *
      tmp.long[["clever.weight"]] *
      tmp.long[["clever.weight.alpha"]] *
      clever *
      ((tmp.long[["delta"]] == process.deltas[jj]) -
        tmp.long[[paste0("at.risk.", name.jj)]] *
          tmp.long[[paste0("P.", name.jj)]] *
          fluct)
    keep <- tmp.long[["at.risk"]] == 1
    data.table(id = tmp.long[["id"]][keep], s = s[keep])[,
      sum(s),
      by = "id"
    ][[2]]
  }

  for (iter in 1:max_iter) {
    if (verbose) {
      message("TMLE iteration ", iter)
    }

    setkey(tmp.long, id, time)

    dt_list <- split(
      tmp.long[,
        !(names(tmp.long) %in%
          c("Q", grep(clever.Q.label, names(tmp.long), value = TRUE))),
        with = FALSE
      ],
      by = "id",
      keep.by = TRUE
    )

    dt_list <- parallel::mclapply(
      dt_list,
      .compute_q_clever,
      states = states,
      process_types = process.types,
      p_prefix = P.prefix,
      parameter = target.name,
      get_years_lost = !is.null(years_lost),
      years_lost_block_size = years_lost,
      clever_by_state = target_by_state,
      mc.cores = .n_cores(cores)
    )

    tmp.long <- rbindlist(dt_list)

    #-- current estimator for target parameter:
    target.est <- mean(tmp.long[, Q[1], by = "id"][[2]])

    if (iter == 1) {
      g.est <- target.est
    }

    eic <- tmp.long[at.risk == 1, Q[1] - target.est, by = "id"][[2]]

    clever.ids <- which(vapply(
      process.names,
      function(process.name) {
        paste0(clever.Q.label, process.name) %in% names(tmp.long)
      },
      logical(1)
    ))

    for (process.jj in clever.ids) {
      eic <- eic + score(process.jj)
    }

    if (iter == 1) {
      se.init <- sqrt(mean(eic^2 / n))
    }

    target.se <- sqrt(mean(eic^2 / n))

    if (one_step) {
      one.step.est <- mean(eic + target.est)
      break
    }

    if (verbose) {
      message("Mean of the efficient influence curve: ", abs(mean(eic)))
    }

    if (iter > min_iter) {
      # se.init (fixed at iteration 1) keeps this threshold stable rather
      # than moving with each iteration's eic.
      if (abs(mean(eic)) <= conv_const * se.init / (log(n))) {
        converged <- TRUE
        if (verbose) {
          message("Converged after ", iter, " iterations.")
        }
        break
      }
    }

    for (process.jj in clever.ids) {
      name.jj <- process.names[process.jj]
      eps.jj <- nleqslv(0, function(eps) mean(score(process.jj, eps)))$x
      if (verbose) {
        message("eps.", name.jj, " = ", eps.jj)
      }
      fluct.cols <- unique(c(
        paste0(P.prefix, name.jj),
        paste0("P.", name.jj)
      ))
      for (col in fluct.cols) {
        if (target_by_state) {
          tmp.long[,
            (col) := get(col) *
              exp(eps.jj * get(paste0(clever.Q.label, name.jj)))
          ]
        } else {
          tmp.long[, (col) := get(col) * exp(eps.jj)]
        }
      }
      for (state.jj in depend.matrix[, unique(state)]) {
        col <- paste0(P.prefix, name.jj, ".", state.jj)
        if (target_by_state) {
          tmp.long[,
            (col) := get(col) *
              exp(eps.jj * get(paste0(clever.Q.label, name.jj, ".", state.jj)))
          ]
        } else {
          tmp.long[, (col) := get(col) * exp(eps.jj)]
        }
      }
    }

    if (length(alpha_list) > 0) {
      .update_alpha_list_weights(tmp.long, alpha_list)
    } else if (z.present) {
      .update_alpha_weights(tmp.long, alpha, z.name)
    }

    names.P <- names(tmp.long)[substr(names(tmp.long), 1, 2) == "P."]

    for (name.P in names.P) {
      if (any(tmp.long[[name.P]] > 1)) {
        if (verbose) {
          message(
            "Transforming ",
            name.P,
            " with 1 - exp(-x) to avoid values > 1."
          )
        }
        tmp.long[, (name.P) := 1 - exp(-tmp.long[[name.P]])]
      }
    }
  }

  if (one_step) {
    out <- list(
      estimate = c(one.step.est = one.step.est, se = target.se, g.est = g.est)
    )
  } else {
    out <- list(
      estimate = c(tmle.est = target.est, se = target.se, g.est = g.est)
    )
  }

  if (output_convergence) {
    out$convergence <- c(
      iter = iter,
      eic.solved.at = abs(mean(eic)),
      converged = converged
    )
  }

  if (length(output_weights) > 0) {
    out$weights <- out.weights
    out$cens.weights <- out.cens.weights
  }

  if (length(output_a_weights) > 0) {
    out$a.weights <- out.a.weights
  }

  if (output_eic) {
    out$eic <- eic
  }

  out
}

# Clever weight of the intervention alpha on process z.name:
# alpha^N_z(t-) * exp(-(alpha - 1) * Lambda_z(t-)).
.update_alpha_weights <- function(tmp.long, alpha, z.name) {
  tmp.long[,
    cum.hazard.z := cumsum(
      get(paste0("at.risk.", z.name)) * get(paste0("P.", z.name))
    ),
    by = "id"
  ]
  tmp.long[, cum.hazard.z.1 := c(0, cum.hazard.z[-.N]), by = "id"]
  tmp.long[,
    clever.weight.alpha := alpha^get(z.name) *
      exp(-(alpha - 1) * cum.hazard.z.1)
  ]
  invisible(tmp.long)
}

# As .update_alpha_weights(), with one time-varying alpha column per process
# in alpha_list.
.update_alpha_list_weights <- function(tmp.long, alpha_list) {
  tmp.long[, clever.weight.alpha := 1]
  for (name.kk in names(alpha_list)) {
    cum.col <- paste0("cum.hazard.", name.kk)
    tmp.long[,
      (cum.col) := cumsum(
        get(paste0("at.risk.", name.kk)) * get(paste0("P.", name.kk))
      ),
      by = "id"
    ]
    tmp.long[, (paste0(cum.col, ".1")) := c(0, get(cum.col)[-.N]), by = "id"]
    tmp.long[,
      clever.weight.alpha := clever.weight.alpha *
        get(paste0("alpha.", name.kk))^get(name.kk) *
        exp(-(get(paste0("alpha.", name.kk)) - 1) * get(paste0(cum.col, ".1")))
    ]
  }
  invisible(tmp.long)
}

.assert_initial_fit <- function(initial_fit) {
  checkmate::assert_list(initial_fit, .var.name = "initial_fit")
  checkmate::assert_names(
    names(initial_fit),
    must.include = c(
      "tmp.long",
      "depend.matrix",
      "process.names",
      "process.deltas",
      "process.types",
      "cens.process.id",
      "at.risks"
    ),
    .var.name = "names(initial_fit)"
  )
}
