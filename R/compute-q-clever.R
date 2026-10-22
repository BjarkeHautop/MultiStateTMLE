# Compute Q and clever covariates for a single id via a discrete-state
# backward recursion.
#
# `dt_id` holds the time-ordered rows of one id, with hazard columns named
# `<p_prefix><name>.<s>` for each process `<name>` and state index `s` in
# 1..S. `states` is the product state space (one column per stateful
# process, plus `state`), used to map state indices and build successor
# mappings. Returns `dt_id` with `Q`, `Q.<s>` and `clever.Q.<name>[0|1]`
# columns added (or `clever.Q.years.lost.<name>` when `get_years_lost`, with
# the integral of the risk taken up to `years_lost_tau`, by default the last
# time of `dt_id`).
.compute_q_clever <- function(
  dt_id,
  states,
  process_types = NULL,
  p_prefix = "P.",
  state_idx_col = "state",
  parameter = "target",
  compute_clever = TRUE,
  get_years_lost = FALSE,
  years_lost_block_size = 10,
  years_lost_tau = NULL,
  clever_by_state = FALSE
) {
  if (is.null(states) || nrow(states) == 0) {
    states <- data.table(state = 1L)
  }
  S <- length(unique(states[[state_idx_col]]))

  Tn <- nrow(dt_id)
  if (Tn < 1) {
    stop("dt_id must have at least one row (time-ordered)")
  }

  # state columns (exclude 'at.risk', the state index and successor columns)
  sstates <- names(states)[substr(names(states), 1, 2) == "s."]
  state_cols_all <- setdiff(
    names(states),
    c("at.risk", state_idx_col, "target.only.in.state", sstates)
  )
  state_processes <- state_cols_all

  time_var_match <- names(states)[substr(names(states), 1, 5) == "time."]

  derived_processes <- state_processes[
    !(state_processes %in% names(process_types))
  ]
  state_processes <- setdiff(state_processes, derived_processes)

  is_target_only_in_state <- "target.only.in.state" %in% names(states)
  if (is_target_only_in_state) {
    target_only_in_state <- states[["target.only.in.state"]]
  }

  # If the declared target is recurrent but not in states, ensure loops still consider it
  if (
    !is.null(process_types) &&
      !is.null(parameter) &&
      parameter %in% names(process_types) &&
      identical(process_types[[parameter]], "recurrent") &&
      !(parameter %in% state_processes)
  ) {
    state_processes <- c(parameter, state_processes)
  }
  # Track which processes actually have columns in 'states'
  stateful_in_states <- intersect(state_processes, state_cols_all)

  # Parse hazard columns in dt_id
  pcols <- grep(
    paste0("^", p_prefix, "[^\\.]+\\.[0-9]+$"),
    names(dt_id),
    value = TRUE
  )
  if (length(pcols) == 0L) {
    warning(
      "No hazard columns found with prefix '",
      p_prefix,
      "'. Returning Q=0"
    )
    out <- copy(dt_id)
    out[, Q := 0]
    return(out[])
  }
  re <- paste0("^", p_prefix, "([^\\.]+)\\.(\\d+)$")
  regs <- regmatches(pcols, regexec(re, pcols))
  found_names <- unique(vapply(
    regs,
    function(x) if (length(x) >= 3) x[2] else NA_character_,
    FUN.VALUE = ""
  ))
  found_names <- found_names[!is.na(found_names)]

  # Build columns_by_name and hazard matrices (Tn x S)
  columns_by_name <- setNames(
    lapply(found_names, function(x) rep(NA_character_, S)),
    found_names
  )
  for (i in seq_along(pcols)) {
    r <- regs[[i]]
    if (length(r) < 3) {
      next
    }
    nm <- r[2]
    sidx <- as.integer(r[3])
    if (is.na(sidx) || sidx < 1L || sidx > S) {
      next
    }
    columns_by_name[[nm]][sidx] <- pcols[i]
  }
  hazard_mats <- list()
  for (nm in found_names) {
    mat <- matrix(0, nrow = Tn, ncol = S)
    for (s in seq_len(S)) {
      colname <- columns_by_name[[nm]][s]
      if (!is.na(colname) && colname %in% names(dt_id)) {
        mat[, s] <- as.numeric(dt_id[[colname]])
      }
    }
    hazard_mats[[nm]] <- mat
  }

  # Partition discovered names into stateful (those in state_processes) and terminal (the rest)
  stateful_names <- intersect(found_names, state_processes)
  terminal_names <- setdiff(found_names, stateful_names)
  # Optionally restrict terminal_names to those declared terminal in process_types
  if (!is.null(process_types)) {
    declared_terminal <- names(process_types)[vapply(
      process_types,
      function(x) identical(x, "terminal"),
      logical(1)
    )]
    if (length(declared_terminal) > 0L) {
      terminal_names <- intersect(terminal_names, declared_terminal)
    }
  }
  stateful_names <- sort(stateful_names)
  terminal_names <- sort(terminal_names)

  # ---- successor mappings
  # Stay mapping: maps row -> state index
  gamma_stay_idx <- unique(states[[state_idx_col]])
  # Jump mapping, aligned with state_processes: read from the s.<process>
  # successor columns of states where present, identity otherwise.
  M <- length(state_processes)
  gamma_jump_idx <- vector("list", M)
  for (j in seq_len(M)) {
    procj_name <- state_processes[j]
    procj_name_s <- grep(procj_name, sstates, value = TRUE)
    if (procj_name %in% stateful_in_states && length(procj_name_s) > 0) {
      gamma_jump_idx[[j]] <- states[[procj_name_s]]
    } else {
      gamma_jump_idx[[j]] <- seq_len(S)
    }
  }

  # ---- classify target
  target_name <- parameter
  if (identical(parameter, "target")) {
    # choose default: first terminal found, else first found name
    if (length(terminal_names) > 0L) {
      target_name <- terminal_names[1L]
    } else {
      target_name <- found_names[1L]
    }
  }
  if (!(target_name %in% c(found_names, state_processes))) {
    stop(
      "parameter must be one of discovered names or a recurrent target declared. Found: ",
      paste(c(found_names, state_processes), collapse = ", ")
    )
  }
  if (!is.null(process_types) && target_name %in% names(process_types)) {
    target_type <- process_types[[target_name]]
  } else {
    target_type <- if (target_name %in% terminal_names) {
      "terminal"
    } else if (target_name %in% stateful_in_states) {
      "recurrent"
    } else {
      "one_jump"
    }
  }
  target_one_jump <- identical(target_type, "one_jump")
  target_terminal <- identical(target_type, "terminal")
  target_recurrent <- identical(target_type, "recurrent")
  target_terminal_or_recurrent <- target_terminal | target_recurrent
  target_name_in_states <- target_name %in% names(states)

  is_recurrent_state_processes <-
    target_recurrent * (state_processes == target_name)

  # ---- backward recursion

  ## Precompute hazard sums (one-time, outside tt-loop)
  if (length(terminal_names) > 0L) {
    Pout_sum_mat <- Reduce(
      `+`,
      lapply(terminal_names, function(nm) hazard_mats[[nm]]),
      init = matrix(0, nrow = Tn, ncol = S)
    )
  } else {
    Pout_sum_mat <- matrix(0, nrow = Tn, ncol = S)
  }
  if (length(stateful_names) > 0L) {
    Pproc_sum_mat <- Reduce(
      `+`,
      lapply(stateful_names, function(nm) hazard_mats[[nm]]),
      init = matrix(0, nrow = Tn, ncol = S)
    )
  } else {
    Pproc_sum_mat <- matrix(0, nrow = Tn, ncol = S)
  }

  hazard_arr <- array(0, dim = c(Tn, S, M))
  for (j in seq_along(state_processes)) {
    nm <- state_processes[j]
    if (nm %in% names(hazard_mats)) {
      hazard_arr[,, j] <- hazard_mats[[nm]]
    }
  }

  gamma_jump_idx_mat <- matrix(NA_integer_, nrow = nrow(states), ncol = M)
  for (j in seq_len(M)) {
    gamma_jump_idx_mat[, j] <- gamma_jump_idx[[j]]
  }

  states1 <- unique(states, by = "state")

  # --- attach rowwise Q
  out <- copy(dt_id)
  s_vec <- as.integer(out[[state_idx_col]])
  s_vec_row_index <- as.integer(out[[paste0(state_idx_col, ".row.index")]])

  all_names <- sort(unique(c(state_processes, terminal_names)))
  no_all <- length(all_names)

  if (length(time_var_match) > 0) {
    # For each distinct combination of time.* values in dt_id, which rows of
    # states are compatible with it (a matching time, or Inf for "not yet").
    state_times <- states[, time_var_match, with = FALSE]
    inf_rows <- lapply(time_var_match, function(time_var) {
      is.infinite(state_times[[time_var]])
    })
    names(inf_rows) <- time_var_match

    uniq_time <- unique(dt_id[, time_var_match, with = FALSE])[, row := 1:.N]
    tt_time <- uniq_time[dt_id, on = time_var_match][["row"]]

    match_mats <- lapply(time_var_match, function(time_var) {
      outer(state_times[[time_var]], uniq_time[[time_var]], "==") |
        inf_rows[[time_var]]
    })
    rows_by_time_mat <- Reduce(`&`, match_mats)

    rows_by_time <- lapply(seq_len(nrow(uniq_time)), function(tm) {
      rows_by_time_mat[, tm]
    })
    names(rows_by_time) <- uniq_time[["row"]]
  }

  compute_q_up_to <- function(Tn) {
    # terminal condition at last row: use hazard if present or zero
    if (target_name %in% names(hazard_mats)) {
      Q_T <- hazard_mats[[target_name]][Tn, ]
      if (is_target_only_in_state) {
        Q_T <- target_only_in_state * Q_T
      }
    } else {
      Q_T <- rep(0, S)
    }
    if (target_one_jump && target_name_in_states) {
      Q_T[as.integer(states1[[target_name]]) >= 1] <- 1
    }

    Q_mat <- matrix(0, nrow = Tn, ncol = S)
    Q_mat[Tn, ] <- Q_T
    Q_row <- numeric(Tn)

    if (Tn > 1) {
      for (tt in (Tn - 1):1) {
        Qn <- Q_mat[tt + 1, ]

        Pstay <- 1 - Pout_sum_mat[tt, ] - Pproc_sum_mat[tt, ]

        Qt <- Pstay * Qn[gamma_stay_idx]

        P_states <- hazard_arr[tt, , ]

        if (length(time_var_match) > 0L) {
          rows_gamma <- rows_by_time[[as.character(tt_time[[tt]])]]
          Qn_jump_mat <- matrix(
            Qn[gamma_jump_idx_mat[rows_gamma, ]],
            nrow = S,
            ncol = M
          )
        } else {
          Qn_jump_mat <- matrix(Qn[gamma_jump_idx_mat], nrow = S, ncol = M)
        }

        Qt <- Qt + rowSums(P_states * Qn_jump_mat)

        # add immediate reward for terminal and recurrent targets
        if (target_terminal_or_recurrent) {
          if (is_target_only_in_state) {
            Qt <- Qt + target_only_in_state * hazard_mats[[target_name]][tt, ]
          } else {
            Qt <- Qt + hazard_mats[[target_name]][tt, ]
          }
        }

        # enforce one-jump terminal states -> Q=1 where already observed
        if (target_one_jump && target_name_in_states) {
          already_idx <- which(as.integer(states1[[target_name]]) >= 1)
          if (length(already_idx) > 0L) Qt[already_idx] <- 1
        }

        Q_mat[tt, ] <- Qt
      }
    }

    for (tt in seq_len(Tn)) {
      Q_row[tt] <- Q_mat[tt, s_vec[tt]]
    }

    if (!target_recurrent && target_name_in_states) {
      target_observed_by_state <- as.integer(states1[[target_name]] >= 1)
    } else {
      target_observed_by_state <- integer(S)
    }

    clever_fun <- function(s_vec) {
      clever0 <- matrix(
        NA_real_,
        nrow = Tn,
        ncol = no_all,
        dimnames = list(NULL, all_names)
      )
      clever1 <- matrix(
        NA_real_,
        nrow = Tn,
        ncol = no_all,
        dimnames = list(NULL, all_names)
      )

      # computed for tt = 1..Tn-1 (using Q_mat[tt+1, ...])
      tt_idx <- seq_len(max(1, Tn - 1))
      m <- (Tn > 1) * length(tt_idx)

      # 1) clever0 rows for tt=1..Tn-1: Q_mat[tt+1, stay state]
      if (m > 0) {
        s_no_vec <- gamma_stay_idx[s_vec[1:m]]
        clever0_vals <- Q_mat[cbind(tt_idx + 1L, s_no_vec)]
        clever0[1:m, ] <- matrix(
          rep(clever0_vals, times = no_all),
          nrow = m,
          ncol = no_all
        )
      }
      # last row (tt = Tn): use Q_mat[Tn, ..]
      clever0[Tn, ] <- Q_mat[Tn, gamma_stay_idx[s_vec[Tn]]]

      # 2) terminal processes:
      if (length(terminal_names) > 0L) {
        term_base <- target_observed_by_state[s_vec]
        if (m > 0) {
          clever1[1:m, terminal_names] <- matrix(
            rep(term_base[1:m], length(terminal_names)),
            nrow = m,
            ncol = length(terminal_names)
          )
        }
        clever1[Tn, terminal_names] <- term_base[Tn]

        # if the target is a terminal and should be set to 1:
        if (target_terminal && target_name %in% terminal_names) {
          if (is_target_only_in_state) {
            clever1[, target_name] <- 1 * target_only_in_state[s_vec[1:(m + 1)]]
          } else {
            clever1[, target_name] <- 1L
          }
        }
      }

      # 3) stateful processes (vectorized)
      if (length(state_processes) > 0L) {
        # Tn x M: successor state indices for the state of each row
        succ_idx_mat <- gamma_jump_idx_mat[s_vec_row_index, , drop = FALSE]

        if (m > 0) {
          # extract Q_mat[tt+1, succ_idx_mat[tt, j]] for all tt, j
          rows_long <- rep(tt_idx + 1L, times = ncol(succ_idx_mat))
          cols_long <- as.integer(as.vector(succ_idx_mat[1:m, , drop = FALSE]))
          clever1_stateful <- matrix(
            Q_mat[cbind(rows_long, cols_long)],
            nrow = m,
            ncol = ncol(succ_idx_mat)
          )
          # add recurrence offsets if present
          if (any(is_recurrent_state_processes != 0L)) {
            clever1_stateful <- clever1_stateful +
              matrix(
                rep(is_recurrent_state_processes, each = m),
                nrow = m,
                ncol = length(is_recurrent_state_processes),
                byrow = FALSE
              )
          }
          clever1[1:m, state_processes] <- clever1_stateful
        }
        # last row: tt = Tn
        clever1[Tn, state_processes] <- Q_mat[
          Tn,
          as.integer(succ_idx_mat[Tn, ])
        ] +
          (is_recurrent_state_processes * 1)
      }

      list(clever1 = clever1, clever0 = clever0)
    }

    if (!compute_clever) {
      return(list(Q_row = Q_row))
    }

    clever_s_obs <- clever_fun(s_vec)
    out_list <- list(
      Q_row = Q_row,
      clever1 = clever_s_obs$clever1,
      clever0 = clever_s_obs$clever0
    )
    for (s in seq_len(S)) {
      if (clever_by_state) {
        clever_s <- clever_fun(rep(s, length = length(s_vec)))
        out_list[[paste0("clever1.", s)]] <- clever_s$clever1
        out_list[[paste0("clever0.", s)]] <- clever_s$clever0
      }
      out_list[[paste0("Q.", s)]] <- Q_mat[, s]
    }
    out_list
  }

  q_out <- compute_q_up_to(Tn)

  out[, Q := q_out$Q_row]

  if (compute_clever) {
    for (nm in all_names) {
      out[[paste0("clever.Q.", nm, "0")]] <- q_out$clever0[, nm]
      out[[paste0("clever.Q.", nm, "1")]] <- q_out$clever1[, nm]
      out[[paste0("clever.Q.", nm)]] <- q_out$clever1[, nm] -
        q_out$clever0[, nm]
      if (clever_by_state) {
        for (s in seq_len(S)) {
          out[[paste0("clever.Q.", nm, ".", s)]] <-
            q_out[[paste0("clever1.", s)]][, nm] -
            q_out[[paste0("clever0.", s)]][, nm]
        }
      }
    }
    for (s in seq_len(S)) {
      out[[paste0("Q.", s)]] <- q_out[[paste0("Q.", s)]]
    }
  }

  if (get_years_lost) {
    # k block size for approximation
    k <- years_lost_block_size
    starts <- unique(c(seq(1L, Tn, by = k), Tn))
    ends <- pmin(starts + k - 1L, Tn)
    rep_idx <- floor((starts + ends) / 2)

    # map every time index to block index 1..length(starts)
    block_of <- findInterval(1:Tn, starts)

    Q_years_lost <- 0
    clever_years_lost <- matrix(
      0,
      nrow = Tn,
      ncol = length(all_names),
      dimnames = list(NULL, all_names)
    )

    if (clever_by_state) {
      clever_years_lost_list <- lapply(seq_len(S), function(s) {
        clever_years_lost
      })
    }

    # The risk is a step function of time, constant from each time of
    # dt_id to the next (or to tau), so its integral weights the risk at
    # each time by the length of the interval that follows it.
    if (is.null(years_lost_tau)) {
      years_lost_tau <- dt_id[["time"]][Tn]
    }
    dt_vec <- diff(c(dt_id[["time"]], years_lost_tau))
    last_block <- NA_integer_
    q_out_block <- NULL

    for (tt in Tn:1) {
      b <- block_of[tt]
      if (is.na(last_block) || b != last_block) {
        # new block encountered -> compute once for block representative
        q_out_block <- compute_q_up_to(rep_idx[b])
        last_block <- b
      }

      # reuse q_out_block for all tt in the same block
      Q_years_lost <- Q_years_lost + dt_vec[tt] * q_out_block$Q_row[1]

      clever_diff <- q_out_block$clever1 - q_out_block$clever0
      imax <- min(tt, nrow(clever_diff))

      clever_years_lost[1:imax, ] <- clever_years_lost[1:imax, ] +
        dt_vec[tt] * out[["at.risk"]][1:imax] * clever_diff[1:imax, ]

      if (clever_by_state) {
        for (s in seq_len(S)) {
          clever_diff <- q_out_block[[paste0("clever1.", s)]] -
            q_out_block[[paste0("clever0.", s)]]
          clever_years_lost_list[[s]][1:imax, ] <-
            clever_years_lost_list[[s]][1:imax, ] +
            dt_vec[tt] * out[["at.risk"]][1:imax] * clever_diff[1:imax, ]
        }
      }
    }

    out[["Q"]] <- Q_years_lost

    for (nm in all_names) {
      out[[paste0("clever.Q.years.lost.", nm)]] <- clever_years_lost[, nm]
      if (clever_by_state) {
        for (s in seq_len(S)) {
          out[[paste0("clever.Q.years.lost.", nm, ".", s)]] <-
            clever_years_lost_list[[s]][, nm]
        }
      }
    }
  }

  out[]
}
