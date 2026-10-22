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
#' @param z_name String. Name of the process in `initial_fit` whose
#'   intensity `alpha` scales. Required unless `alpha = 1` or `alpha_list`
#'   is given.
#' @param alpha_list Experimental, and may change. Optional named list of
#'   functions, one per process, giving a time-varying multiplicative
#'   intervention on that process's intensity. Each function is called with
#'   the columns of `initial_fit$tmp_long` named by its arguments, e.g.
#'   `time`, a baseline covariate such as `L0`, or a process such as `z`
#'   (its number of jumps so far), and must return one non-negative factor
#'   per row. For example, `list(z = function(time) ifelse(time < 0.5, 0.5,
#'   1))` halves the intensity of `z` before time 0.5. Supersedes `alpha`
#'   and `z_name`.
#' @param a Optional baseline treatment value (`0` or `1`) to intervene on.
#'   Must be one of the values prepared by [fit_initial()], `initial_fit$a`.
#' @param estimand The target parameter: `"risk"`, the risk of `target` by
#'   `tau` (or the expected number of events for a recurrent `target`), or
#'   `"years_lost"`, the expected time lost, i.e. the integral of the risk
#'   of `target` up to `tau`.
#' @param block_size Positive integer. For `estimand = "years_lost"`, the
#'   integral is approximated by evaluating the risk once per block of
#'   `block_size` consecutive time points. `1` gives the exact integral of
#'   the estimated risk; the computation time is roughly proportional to
#'   `1 / block_size`.
#' @param only_first Logical. If `TRUE`, only the first `target` event
#'   counts, and only the first jump of `z_name` is intervened on.
#' @param target_only_in_state Experimental, and may change. Optional
#'   function of the state space `initial_fit$depend_matrix`, a data frame
#'   with one row per state and a column per process that other processes
#'   depend on (its number of jumps), returning a logical vector marking the
#'   states in which `target` events count. For example,
#'   `function(states) states[["z"]] == 0` counts only deaths before
#'   illness.
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
#' @param cores Positive integer. Number of cores used by
#'   [parallel::mclapply()] for the per-individual computations.
#' @param verbose Logical. If `TRUE`, print progress.
#'
#' @return An object of class `tmle_alpha`, a list with elements:
#'   * `estimate`: a `data.table` with columns `target`, `estimate`, its
#'     standard error `se`, and the bounds `lower` and `upper` of a 95% Wald
#'     confidence interval.
#'   * `initial_estimate`: the initial plug-in estimate, before targeting.
#'   * `method`: `"tmle"`, or `"one_step"` if `one_step = TRUE`.
#'   * `estimand`: as given.
#'   * `alpha`, `z_name` and `a`: the intervention (`alpha` is `NULL` if
#'     `alpha_list` is given).
#'   * `eic`: a list, named by `target`, with the efficient influence curve,
#'     one value per individual.
#'   * `convergence`: a named vector with the number of iterations `iter`,
#'     the absolute mean of the efficient influence curve `eic_solved_at`,
#'     and whether the iterations `converged`. A warning is given if they
#'     did not.
#'   * `weights`: a matrix with the quantiles (columns `q50`, `q90`, `q99`
#'     and `q100`) over individuals of each individual's maximal
#'     `intervention`, `censoring` and `treatment` weight, and of their
#'     product `total`, the weight entering the targeting step (rows), as a
#'     positivity diagnostic. The intervention weights are those after
#'     truncation, and the treatment weights are 1 when `a` is not given.
#'   * `n_truncated`: the number of individuals whose initial intervention
#'     weight was truncated (0 if `truncate_weights = 0`).
#'
#' @examples
#' fit <- fit_initial(
#'   illness_death,
#'   fit_types = list(
#'     z = list(
#'       delta = 2,
#'       model = ~ L0,
#'       at_risk = function(dt) dt[["z"]] == 0
#'     ),
#'     death = list(delta = 1, model = ~ L0 + z),
#'     censoring = list(delta = 0, model = ~ L0)
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
  z_name = NULL,
  alpha_list = NULL,
  a = NULL,
  estimand = c("risk", "years_lost"),
  block_size = 10,
  only_first = FALSE,
  target_only_in_state = NULL,
  target_by_state = FALSE,
  conv_const = 1,
  one_step = FALSE,
  max_iter = 10,
  min_iter = 0,
  truncate_weights = 0,
  cores = 1,
  verbose = FALSE
) {
  .assert_initial_fit(initial_fit)
  .assert_a(a, initial_fit)
  checkmate::assert_string(target)
  checkmate::assert_number(alpha, lower = 0, finite = TRUE)
  checkmate::assert_string(z_name, null.ok = TRUE)
  checkmate::assert_list(
    alpha_list,
    types = "function",
    names = "unique",
    null.ok = TRUE
  )
  estimand <- match.arg(estimand)
  checkmate::assert_count(block_size, positive = TRUE)
  checkmate::assert_flag(only_first)
  checkmate::assert_function(target_only_in_state, null.ok = TRUE)
  checkmate::assert_flag(target_by_state)
  checkmate::assert_number(conv_const, lower = 0, finite = TRUE)
  checkmate::assert_flag(one_step)
  checkmate::assert_count(max_iter, positive = TRUE)
  checkmate::assert_count(min_iter)
  checkmate::assert_number(truncate_weights, lower = 0, finite = TRUE)
  checkmate::assert_count(cores, positive = TRUE)
  checkmate::assert_flag(verbose)

  tmp.long <- copy(initial_fit$tmp_long)
  depend.matrix <- copy(initial_fit$depend_matrix)

  n <- length(unique(tmp.long[["id"]]))

  process.names <- initial_fit$process_names
  process.deltas <- initial_fit$process_deltas
  process.types <- initial_fit$process_types

  at.risks <- initial_fit$at_risks

  target.id <- which(tolower(process.names) == tolower(target))
  if (length(target.id) != 1) {
    stop(
      "target must be one of the processes in initial_fit: ",
      paste(process.names, collapse = ", "),
      "."
    )
  }
  target.name <- process.names[target.id]

  if (is.null(z_name)) {
    if (length(alpha_list) == 0 && !identical(alpha, 1)) {
      stop("z_name is required when alpha != 1.")
    }
    z.process.id <- integer(0)
  } else {
    z.process.id <- which(tolower(process.names) == tolower(z_name))
    if (length(z.process.id) != 1) {
      stop(
        "z_name must be one of the processes in initial_fit: ",
        paste(process.names, collapse = ", "),
        "."
      )
    }
  }

  if (length(z.process.id) == 1) {
    z.name <- process.names[z.process.id]
  } else {
    z.process.id <- Inf
    z.name <- "not.applied"
  }
  z.present <- length(alpha_list) > 0 || is.finite(z.process.id)

  #--------------------------------
  #-- initialization step:

  states <- copy(depend.matrix)
  state.names <- names(states)
  setkeyv(states, c("state", state.names[state.names %in% process.names]))

  if (only_first) {
    process.types[[target.name]] <- "one_jump"
    tmp.long <- tmp.long[tmp.long[[target.name]] <= 1]
  }

  tmp.long[, at.risk := (time <= final.time)]

  if (target.name %in% names(at.risks)) {
    if (process.types[target.name] == "one_jump") {
      tmp.long[, at.risk := at.risk * (at.risks[[target.name]](tmp.long))]
    }
    states[, at.risk := at.risks[[target.name]](depend.matrix)]
  }

  if (only_first && is.finite(z.process.id)) {
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

  # Recompute the intervention's clever weights from the current hazards,
  # truncated at truncate_weights. Returns the number of individuals whose
  # weight was truncated.
  update_weights <- function() {
    if (length(alpha_list) > 0) {
      .update_alpha_list_weights(tmp.long, alpha_list)
    } else {
      .update_alpha_weights(tmp.long, alpha, z.name)
    }
    if (truncate_weights == 0) {
      return(0)
    }
    truncated <- tmp.long[["time"]] <= tmp.long[["final.time"]] &
      tmp.long[["clever.weight.alpha"]] > truncate_weights
    tmp.long[truncated, clever.weight.alpha := truncate_weights]
    length(unique(tmp.long[["id"]][truncated]))
  }

  if (z.present) {
    n.truncated <- update_weights()
  } else {
    tmp.long[, clever.weight.alpha := 1]
    n.truncated <- 0
  }

  # Quantiles over individuals of each individual's maximal weight up to
  # their final time, as a positivity diagnostic.
  weight_quantiles <- function(w) {
    keep <- tmp.long[["time"]] <= tmp.long[["final.time"]]
    per_id <- tapply(w[keep], tmp.long[["id"]][keep], max)
    quantile(per_id, c(0.5, 0.9, 0.99, 1), names = FALSE)
  }

  # Taken before clever.weight is overwritten by the combined censoring and
  # treatment weight below.
  cens.weights <- weight_quantiles(tmp.long[["clever.weight"]])

  if (length(a) > 0) {
    a.weight.col <- paste0("clever.weight.a", a)
    a.weights <- weight_quantiles(tmp.long[,
      (A0.obs == a) / (pi.A0.1^(A0.obs == 1) * (1 - pi.A0.1)^(A0.obs == 0))
    ])
    tmp.long[, clever.weight := tmp.long[[a.weight.col]]]
  } else {
    a.weights <- rep(1, 4)
  }

  weights <- rbind(
    intervention = weight_quantiles(tmp.long[["clever.weight.alpha"]]),
    censoring = cens.weights,
    treatment = a.weights,
    total = weight_quantiles(
      tmp.long[["clever.weight"]] * tmp.long[["clever.weight.alpha"]]
    )
  )
  colnames(weights) <- c("q50", "q90", "q99", "q100")

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
  years_lost <- estimand == "years_lost"
  clever.Q.label <- if (years_lost) {
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
      get_years_lost = years_lost,
      years_lost_block_size = block_size,
      years_lost_tau = initial_fit$tau,
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

    if (z.present) {
      update_weights()
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

  if (!one_step && !converged) {
    warning(
      "TMLE did not converge in max_iter = ",
      max_iter,
      " iterations; the mean of the efficient influence curve is ",
      signif(abs(mean(eic)), 3),
      "."
    )
  }

  structure(
    list(
      estimate = .wald_table(
        target.name,
        if (one_step) one.step.est else target.est,
        target.se
      ),
      initial_estimate = g.est,
      method = if (one_step) "one_step" else "tmle",
      estimand = estimand,
      alpha = if (length(alpha_list) > 0) NULL else alpha,
      z_name = if (is.finite(z.process.id)) z.name else NULL,
      a = a,
      eic = stats::setNames(list(eic), target.name),
      convergence = c(
        iter = iter,
        eic_solved_at = abs(mean(eic)),
        converged = converged
      ),
      weights = weights,
      n_truncated = n.truncated
    ),
    class = "tmle_alpha"
  )
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
  checkmate::assert_class(initial_fit, "initial_fit")
}

# A baseline treatment value `a` must be one that fit_initial() prepared.
.assert_a <- function(a, initial_fit, var_name = "a") {
  checkmate::assert_int(
    a,
    lower = 0,
    upper = 1,
    null.ok = TRUE,
    .var.name = var_name
  )
  if (is.null(a) || a %in% initial_fit$a) {
    return(invisible(a))
  }
  if (length(initial_fit$a) == 0) {
    stop(
      var_name,
      " requires a treatment_model in fit_initial().",
      call. = FALSE
    )
  }
  stop(
    var_name,
    " must be one of the values prepared by fit_initial(): ",
    paste(initial_fit$a, collapse = ", "),
    ".",
    call. = FALSE
  )
}
