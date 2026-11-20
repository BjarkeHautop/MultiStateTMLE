#' Targeted estimation under a stochastic intensity intervention
#'
#' Estimates a target parameter, such as the risk of death by the horizon
#' `tau` of the initial fit, under a stochastic intensity intervention,
#' optionally combined with an intervention `a` on the baseline treatment.
#' The initial fit from [fit_initial()] is updated by iterative targeted
#' maximum likelihood estimation (TMLE).
#'
#' The intervention is either a fixed intensity scaling,
#' [scale_intensity()], or a calibrated one, [calibrate()], whose scaling
#' is chosen to bring an exposure to a given value. Under a calibrated
#' intervention the target is a composite parameter, and its standard error
#' accounts for the scaling being estimated.
#'
#' @param initial_fit An initial fit, as returned by [fit_initial()].
#' @param target The target parameter: [risk()] or [years_lost()] of one or
#'   more processes, each estimated separately under the same intervention.
#' @param intervention The intervention: [scale_intensity()],
#'   [calibrate()], or `NULL` (default) for none.
#' @param a Optional baseline treatment value (`0` or `1`) to intervene on.
#'   Must be one of the `arms` prepared by [fit_initial()]. Without it, each
#'   individual keeps their observed treatment.
#' @param control A named list of settings for how the estimate is
#'   computed:
#'   * `one_step` (`FALSE`): if `TRUE`, return the one-step estimator
#'     instead of iterating the TMLE update.
#'   * `block_size` (`10`): positive integer. For [years_lost()], the
#'     integral is approximated by evaluating the risk once per block of
#'     `block_size` consecutive time points. `1` gives the exact integral of
#'     the estimated risk; the computation time is roughly proportional to
#'     `1 / block_size`.
#'   * `target_by_state` (`FALSE`): if `TRUE`, use state-specific clever
#'     covariates in the targeting step.
#'   * `conv_const` (`1`): positive number. The TMLE iterations stop when the
#'     absolute mean of the efficient influence curve falls below
#'     `conv_const * se / log(n)`.
#'   * `max_iter` (`10`): maximum number of TMLE iterations.
#'   * `min_iter` (`0`): minimum number of TMLE iterations before checking
#'     convergence.
#'   * `truncate_weights` (`Inf`): positive number at which the clever
#'     weights of the intervention are truncated. `Inf` is no truncation.
#' @param cores Positive integer. Number of cores used by
#'   [parallel::mclapply()] for the per-individual computations.
#' @param verbose Logical. If `TRUE`, print progress.
#'
#' @return An object of class `tmle_estimate`, a list with elements:
#'   * `estimate`: a `data.table` with one row per target process and
#'     columns `target`, `estimate`, its standard error `se`, and the bounds
#'     `lower` and `upper` of a 95% Wald confidence interval. Under a
#'     calibrated intervention, `se_crude` is the standard error that treats
#'     the scaling as known.
#'   * `target`, `intervention` and `a`: as given.
#'   * `method`: `"tmle"`, or `"one_step"` if `control$one_step = TRUE`.
#'   * `eic`: a list, named by target process, with the efficient influence
#'     curve, one value per individual.
#'   * `checks`: a list of diagnostics:
#'     * `iter`, `eic_solved_at` and `converged`: vectors, named by target
#'       process, of the number of iterations, the absolute mean of the
#'       efficient influence curve, and whether the iterations converged. A
#'       warning is given if they did not. Each iteration checks convergence
#'       before updating, so a converged fit made `iter - 1` targeting steps:
#'       `iter = 1` means the initial fit already met the criterion and was
#'       not updated.
#'     * `weights`: a matrix with the quantiles (columns `q50`, `q90`, `q99`
#'       and `q100`) over individuals of each individual's maximal
#'       `intervention`, `censoring` and `treatment` weight, and of their
#'       product `total`, the weight entering the targeting step (rows), as
#'       a positivity diagnostic. The intervention weights are those after
#'       truncation, and the treatment weights are 1 when `a` is not given.
#'     * `n_truncated`: the number of individuals whose initial intervention
#'       weight was truncated (0 if `control$truncate_weights = Inf`).
#'
#'     `weights` and `n_truncated` are those of the first target process.
#'
#'   Without a calibrated intervention it also has `initial_estimate`, the
#'   initial plug-in estimates before targeting, named by target process.
#'   With one, it has:
#'   * `calibration`: a `data.table` with the estimated scaling `alpha` and
#'     the exposure `theta` it attains, with standard errors and confidence
#'     intervals.
#'   * `eic_alpha_known`: as `eic`, treating the scaling as known.
#'   * `alpha_eic`: the efficient influence curve of the scaling.
#'   * `checks$calibration`: a list with `converged` (whether a scaling
#'     attaining `theta` was found), `dist` (how far from `theta` it is),
#'     `c_n` (the tolerance), and the estimated derivatives with respect to
#'     the scaling of the exposure, `derivative_exposure`, and of the
#'     targets, `derivative_target`.
#'
#' @examples
#' fit <- fit_initial(
#'   illness_death,
#'   models = list(
#'     hazard_model(illness ~ biomarker + treatment, limit = 1),
#'     death ~ biomarker + treatment + illness,
#'     censoring ~ biomarker
#'   ),
#'   tau = 1
#' )
#'
#' # The risk of death by tau = 1 without intervention, and if the intensity
#' # of illness were halved:
#' observed <- tmle_estimate(fit, risk(death))
#' halved <- tmle_estimate(fit, risk(death), scale_intensity(illness = 0.5))
#' observed
#' halved
#'
#' # Several targets at once, here the risks of death and of illness themselves:
#' both <- tmle_estimate(
#'   fit,
#'   risk(death, illness),
#'   scale_intensity(illness = 0.5)
#' )
#' both$estimate
#'
#' # The expected years lost before tau rather than the risk:
#' coef(tmle_estimate(fit, years_lost(death), scale_intensity(illness = 0.5)))
#'
#' # Diagnostics: the quantiles over individuals of the largest weights, as a
#' # check of positivity.
#' halved$checks$weights
#'
#' \donttest{
#' # A calibrated intervention: the risk of death if the risk of illness fell
#' # to 80% of its observed value. The first table gives the scaling alpha that
#' # achieves it and the risk of illness, theta, it attains:
#' cal <- tmle_estimate(
#'   fit,
#'   risk(death),
#'   calibrate(risk(illness), relative_change(0.8))
#' )
#' cal
#'
#' # Its standard error accounts for alpha being estimated; se_crude treats
#' # alpha as known:
#' cal$estimate[, c("target", "estimate", "se", "se_crude")]
#' }
#' @export
tmle_estimate <- function(
  initial_fit,
  target,
  intervention = NULL,
  a = NULL,
  control = list(),
  cores = 1,
  verbose = FALSE
) {
  .assert_initial_fit(initial_fit)
  .assert_target(target, initial_fit)
  .assert_intervention(intervention, initial_fit)
  .assert_a(a, initial_fit)
  control <- .tmle_control(control)
  checkmate::assert_count(cores, positive = TRUE)
  checkmate::assert_flag(verbose)

  if (!is.null(intervention) && intervention$type == "calibrate") {
    return(.tmle_calibrated(
      initial_fit,
      target,
      intervention,
      a,
      control,
      cores,
      verbose
    ))
  }

  fits <- lapply(target$processes, function(process) {
    .tmle_fixed(
      initial_fit,
      process,
      target,
      intervention,
      a,
      control,
      cores,
      verbose
    )
  })
  out <- .combine_targets(fits)
  out$target <- target
  out$intervention <- intervention
  out$a <- a
  structure(
    out[c(
      "estimate",
      "initial_estimate",
      "target",
      "intervention",
      "a",
      "method",
      "eic",
      "checks"
    )],
    class = "tmle_estimate"
  )
}

# The TMLE of `target` for the single target process `process` under the
# fixed intervention `intervention` (scale_intensity() or NULL). Returns a
# list with elements estimate, initial_estimate, method, eic and checks.
.tmle_fixed <- function(
  initial_fit,
  process,
  target,
  intervention,
  a,
  control,
  cores,
  verbose
) {
  alpha <- if (is.null(intervention)) list() else intervention$alpha

  one_step <- control$one_step
  block_size <- control$block_size
  target_by_state <- control$target_by_state
  conv_const <- control$conv_const
  max_iter <- control$max_iter
  min_iter <- control$min_iter
  truncate_weights <- control$truncate_weights

  tmp.long <- copy(initial_fit$tmp_long)
  depend.matrix <- copy(initial_fit$depend_matrix)
  .set_treatment(tmp.long, initial_fit, a)

  n <- length(unique(tmp.long[["id"]]))

  process.names <- initial_fit$process_names
  process.deltas <- initial_fit$process_deltas
  process.types <- initial_fit$process_types

  at.risks <- initial_fit$at_risks

  target.name <- process
  z.present <- length(alpha) > 0

  #--------------------------------
  #-- initialization step:

  states <- copy(depend.matrix)
  state.names <- names(states)
  setkeyv(states, c("state", state.names[state.names %in% process.names]))

  # Only the first event of a recurrent target counts: as a one-jump target,
  # its Q is 1 in every state after its first event, so the rows there add
  # nothing to the clever covariates. They are kept, since the Q recursion
  # of each individual runs over all their rows up to tau.
  if (target$only_first && process.types[[target.name]] == "recurrent") {
    if (!target.name %in% names(tmp.long)) {
      stop(
        "only_first needs the number of events of ",
        target.name,
        " in the state space; use it as a covariate in a hazard model.",
        call. = FALSE
      )
    }
    process.types[[target.name]] <- "one_jump"
  }

  tmp.long[, at.risk := (time <= final.time)]

  if (target.name %in% names(at.risks)) {
    if (process.types[target.name] == "one_jump") {
      tmp.long[, at.risk := at.risk * (at.risks[[target.name]](tmp.long))]
    }
    states[, at.risk := at.risks[[target.name]](depend.matrix)]
  }

  if (!is.null(target$when)) {
    when <- .state_function(
      target$when[[2]],
      setdiff(names(states), "state"),
      environment(target$when),
      "`when`"
    )
    states[, target.only.in.state := when(states)]
  }

  for (process.name in process.names[!(process.names %in% names(at.risks))]) {
    tmp.long[, (paste0("at.risk.", process.name)) := at.risk]
  }

  P.prefix <- if (length(a) > 0) paste0("P.a", a, ".") else "P."

  for (name.kk in names(alpha)) {
    if (
      !name.kk %in% names(tmp.long) &&
        process.types[[name.kk]] != "terminal"
    ) {
      stop(
        "To intervene on ",
        name.kk,
        ", its number of events must be in the state space; use it as a ",
        "covariate in a hazard model.",
        call. = FALSE
      )
    }
  }

  # The factor at the observed counts, for the weights and scores, with the
  # counts capped as in the state space so that it matches the factor of
  # the observed state below.
  count.names <- intersect(process.names, names(depend.matrix))
  data.observed <- as.list(tmp.long)
  for (count.name in count.names) {
    data.observed[[count.name]] <- pmin(
      data.observed[[count.name]],
      max(depend.matrix[[count.name]])
    )
  }
  for (name.kk in names(alpha)) {
    tmp.long[,
      (paste0("alpha.", name.kk)) := .alpha_values(
        alpha[[name.kk]],
        data.observed,
        name.kk
      )
    ]
  }

  #--------------------------------
  #-- compute clever weights:

  # Recompute the intervention's clever weights from the current hazards,
  # truncated at truncate_weights. Returns the number of individuals whose
  # weight was truncated.
  update_weights <- function() {
    .update_alpha_weights(tmp.long, names(alpha))
    if (is.infinite(truncate_weights)) {
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

  # The hazard of state s is scaled by the factor evaluated at the counts of
  # state s, not at the observed ones, since the Q recursion moves through
  # all states.
  for (name.kk in names(alpha)) {
    for (state.kk in depend.matrix[, unique(state)]) {
      varname <- paste0(P.prefix, name.kk, ".", state.kk)
      if (!varname %in% names(tmp.long)) {
        next
      }
      counts <- depend.matrix[state == state.kk][1]
      data.kk <- as.list(tmp.long)
      for (count.name in count.names) {
        data.kk[[count.name]] <- rep(counts[[count.name]], nrow(tmp.long))
      }
      tmp.long[,
        (varname) := .alpha_values(
          alpha[[name.kk]],
          data.kk,
          name.kk
        ) *
          get(varname)
      ]
    }
  }

  #--------------------------------
  #-- tmle iterations:

  converged <- FALSE
  years_lost <- target$estimand == "years_lost"
  clever.Q.label <- if (years_lost) {
    "clever.Q.years.lost."
  } else {
    "clever.Q."
  }

  # The intervention's factor on the score of process number jj: its alpha
  # for the intervened processes, else 1.
  alpha.factor <- function(jj) {
    name.jj <- process.names[jj]
    if (name.jj %in% names(alpha)) {
      return(tmp.long[[paste0("alpha.", name.jj)]])
    }
    1
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
          message("Converged after ", iter - 1, " targeting step(s).")
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

    names.P <- names(tmp.long)[startsWith(names(tmp.long), "P.")]

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
      "TMLE for target '",
      target.name,
      "' did not converge in control$max_iter = ",
      max_iter,
      " iterations; the mean of the efficient influence curve is ",
      signif(abs(mean(eic)), 3),
      "."
    )
  }

  list(
    estimate = .wald_table(
      target.name,
      if (one_step) one.step.est else target.est,
      target.se
    ),
    initial_estimate = stats::setNames(g.est, target.name),
    method = if (one_step) "one_step" else "tmle",
    eic = stats::setNames(list(eic), target.name),
    checks = list(
      iter = stats::setNames(iter, target.name),
      eic_solved_at = stats::setNames(abs(mean(eic)), target.name),
      converged = stats::setNames(converged, target.name),
      weights = weights,
      n_truncated = n.truncated
    )
  )
}

# Combine single-process .tmle_fixed() fits under the same intervention into
# one. The weight diagnostics are those of the first process.
.combine_targets <- function(fits) {
  out <- fits[[1]]
  out$estimate <- rbindlist(lapply(fits, function(fit) fit$estimate))
  out$initial_estimate <- unlist(lapply(fits, function(fit) {
    fit$initial_estimate
  }))
  out$eic <- unlist(lapply(fits, function(fit) fit$eic), recursive = FALSE)
  for (check in c("iter", "eic_solved_at", "converged")) {
    out$checks[[check]] <- unlist(lapply(fits, function(fit) {
      fit$checks[[check]]
    }))
  }
  out
}

# The factor `alpha` of scale_intensity() on process `name`, one value per
# row of `data` (a data.table or a list of equal-length columns): a number,
# or a one-sided formula evaluated in the columns of `data`.
.alpha_values <- function(alpha, data, name) {
  n.rows <- length(data[["time"]])
  if (!inherits(alpha, "formula")) {
    return(rep(alpha, n.rows))
  }
  values <- tryCatch(
    eval(alpha[[2]], data, environment(alpha)),
    error = function(e) {
      stop(
        "Could not evaluate the factor of ",
        name,
        " in scale_intensity(), ",
        .deparse_one_line(alpha),
        ": ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  if (length(values) == 1) {
    values <- rep(values, n.rows)
  }
  if (
    !is.numeric(values) ||
      length(values) != n.rows ||
      anyNA(values) ||
      any(values < 0)
  ) {
    stop(
      "The factor of ",
      name,
      " in scale_intensity() must give one non-negative number per row.",
      call. = FALSE
    )
  }
  values
}

# Clever weight of the intervention on the processes `names`, each with its
# factor in the column alpha.<name>: the likelihood ratio of the intervened
# to the observed process up to t-, the product over them and over earlier
# rows of alpha^dN * exp(-(alpha - 1) * dLambda). dN is read off the count
# column of the process; a terminal process without one has no event while
# follow-up continues.
.update_alpha_weights <- function(tmp.long, names) {
  tmp.long[, clever.weight.alpha := 1]
  for (name.kk in names) {
    alpha.kk <- tmp.long[[paste0("alpha.", name.kk)]]
    d.lambda <- tmp.long[[paste0("at.risk.", name.kk)]] *
      tmp.long[[paste0("P.", name.kk)]]
    if (name.kk %in% names(tmp.long)) {
      tmp.long[, .d.n := c(diff(get(name.kk)), 0), by = "id"]
      d.n <- tmp.long[[".d.n"]]
      tmp.long[, .d.n := NULL]
    } else {
      d.n <- 0
    }
    tmp.long[, .ratio := alpha.kk^d.n * exp(-(alpha.kk - 1) * d.lambda)]
    tmp.long[,
      clever.weight.alpha := clever.weight.alpha * c(1, cumprod(.ratio)[-.N]),
      by = "id"
    ]
    tmp.long[, .ratio := NULL]
  }
  invisible(tmp.long)
}

# A scale_intensity() of process `name` by the number `alpha`.
.scale_one <- function(name, alpha) {
  do.call(scale_intensity, stats::setNames(list(alpha), name))
}

.assert_initial_fit <- function(initial_fit) {
  checkmate::assert_class(initial_fit, "initial_fit")
}

# With a treatment_model, fit_initial() keeps the observed treatment in
# A0.obs, predicts the hazards under each arm `aa` as P.a<aa>.* columns, and
# leaves the treatment column at the last arm. Set the treatment column to
# `a`, or, without `a`, to the observed treatment, with each individual's
# hazards (P.*) those of their observed arm, so that the estimate is under
# the observed treatment as for a fit without a treatment_model.
.set_treatment <- function(tmp.long, initial_fit, a) {
  if (is.null(initial_fit$treatment_model)) {
    return(invisible(tmp.long))
  }
  treatment.name <- all.vars(initial_fit$treatment_model[[2]])
  if (!is.null(a)) {
    set(tmp.long, j = treatment.name, value = a)
    return(invisible(tmp.long))
  }
  if (!setequal(initial_fit$arms, 0:1)) {
    stop(
      "Without a, tmle_estimate() estimates under the observed treatment, ",
      "which needs fit_initial() to prepare both arms, but it prepared only ",
      "a = ",
      paste(initial_fit$arms, collapse = ", "),
      ". Give a, or refit with arms = c(0, 1).",
      call. = FALSE
    )
  }
  treated <- tmp.long[["A0.obs"]] == 1
  for (col1 in grep("^P[.]a1[.]", names(tmp.long), value = TRUE)) {
    rest <- substring(col1, nchar("P.a1.") + 1)
    col0 <- paste0("P.a0.", rest)
    set(
      tmp.long,
      j = paste0("P.", rest),
      value = ifelse(treated, tmp.long[[col1]], tmp.long[[col0]])
    )
  }
  set(tmp.long, j = treatment.name, value = tmp.long[["A0.obs"]])
  invisible(tmp.long)
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
  if (is.null(a) || a %in% initial_fit$arms) {
    return(invisible(a))
  }
  if (length(initial_fit$arms) == 0) {
    stop(
      var_name,
      " requires a treatment_model in fit_initial().",
      call. = FALSE
    )
  }
  stop(
    var_name,
    " must be one of the values prepared by fit_initial(): ",
    paste(initial_fit$arms, collapse = ", "),
    ".",
    call. = FALSE
  )
}
