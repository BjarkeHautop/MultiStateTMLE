#' Solve for the intervention achieving a target value
#'
#' Finds the intensity intervention `alpha` at which the estimate returned
#' by `fun(alpha)` equals `theta`, assuming the estimate is increasing in
#' `alpha`. The search first expands a bracket around the solution, then
#' narrows it by secant steps (with bisection as a fallback), caching
#' evaluations of `fun` along the way.
#'
#' This solves the calibration equation of Rytgaard and van der Laan, and is
#' used by [calibrate_alpha()].
#'
#' @param fun A function of `alpha` returning a [tmle_alpha()] fit, or a
#'   list like it with elements `estimate` (a data frame with columns
#'   `estimate` and `se`) and `eic` (a list with the influence curve).
#' @param theta Number. The target value of the estimate.
#' @param c_n Positive number. Tolerance; the search stops once the estimate
#'   is within `c_n` of `theta`.
#' @param alpha_init Number. Starting value of `alpha`.
#' @param expand_up Number greater than 1. Factor by which `alpha` is
#'   multiplied while searching upwards for a bracket.
#' @param expand_down Number between 0 and 1. Factor by which `alpha` is
#'   multiplied while searching downwards for a bracket.
#' @param max_iter Positive integer. Maximum number of steps in each of the
#'   bracketing and narrowing phases.
#' @param alpha_min,alpha_max Numbers with `0 < alpha_min < alpha_max`.
#'   Bounds on `alpha`.
#' @param verbose Logical. If `TRUE`, print progress every `trace_every`
#'   evaluations.
#' @param trace_every Positive integer. How often to print progress when
#'   `verbose = TRUE`.
#'
#' @return A list with elements `alpha_hat` (the best `alpha` found),
#'   `converged` (whether it is within `c_n` of `theta`), `dist` (its
#'   distance to `theta`), `c_n`, `grid` (a data frame of all evaluated
#'   `alpha`, `psi` and `se`) and `eic` (the influence curve at `alpha_hat`).
#'
#' @examples
#' # psi(alpha) = 0.6 * (1 - exp(-alpha)) has a known inverse, so
#' # solve_alpha() should recover alpha = -log(1 - 0.3 / 0.6).
#' fun <- function(alpha) {
#'   psi <- 0.6 * (1 - exp(-alpha))
#'   list(estimate = data.frame(estimate = psi, se = 0.01), eic = list(0))
#' }
#' res <- solve_alpha(fun, theta = 0.3, c_n = 1e-6)
#' res$alpha_hat
#' -log(1 - 0.3 / 0.6)
#' @export
solve_alpha <- function(
  fun,
  theta,
  c_n,
  alpha_init = 1,
  expand_up = 2,
  expand_down = 0.25,
  max_iter = 100,
  alpha_min = 1e-3,
  alpha_max = 100,
  verbose = FALSE,
  trace_every = 1L
) {
  checkmate::assert_function(fun)
  checkmate::assert_number(theta, finite = TRUE)
  checkmate::assert_number(c_n, lower = 0, finite = TRUE)
  checkmate::assert_number(alpha_init, finite = TRUE)
  checkmate::assert_number(expand_up, lower = 1, finite = TRUE)
  checkmate::assert_number(expand_down, lower = 0, upper = 1)
  checkmate::assert_count(max_iter, positive = TRUE)
  checkmate::assert_number(alpha_min, lower = 0, finite = TRUE)
  checkmate::assert_number(alpha_max, finite = TRUE)
  if (alpha_min <= 0 || alpha_max <= alpha_min) {
    stop("Require 0 < alpha_min < alpha_max.")
  }
  checkmate::assert_flag(verbose)
  checkmate::assert_count(trace_every, positive = TRUE)

  # A little extra room for expansion + bracketing
  max_evals <- 2L * max_iter + 5L

  alpha_hist <- numeric(max_evals)
  psi_hist <- numeric(max_evals)
  se_hist <- numeric(max_evals)
  eic_hist <- vector("list", max_evals)

  n_eval <- 0L
  best_idx <- NA_integer_
  best_dist <- Inf

  cache <- new.env(parent = emptyenv())

  record_eval <- function(rec) {
    n_eval <<- n_eval + 1L
    alpha_hist[n_eval] <<- rec$alpha
    psi_hist[n_eval] <<- rec$psi
    se_hist[n_eval] <<- rec$se
    eic_hist[[n_eval]] <<- rec$eic

    d <- abs(rec$psi - theta)
    if (d < best_dist) {
      best_dist <<- d
      best_idx <<- n_eval
    }
    invisible(NULL)
  }

  eval_psi <- function(alpha) {
    alpha <- max(min(alpha, alpha_max), alpha_min)
    key <- sprintf("%.15g", alpha)

    if (exists(key, envir = cache, inherits = FALSE)) {
      out <- get(key, envir = cache, inherits = FALSE)
      out$cached <- TRUE
      return(out)
    }

    psi <- fun(alpha = alpha)

    out <- list(
      alpha = alpha,
      psi = .point_estimate(psi),
      se = as.numeric(psi[["estimate"]][["se"]][1]),
      eic = .eic(psi),
      cached = FALSE
    )

    assign(key, out, envir = cache)
    out
  }

  progress <- function(
    stage,
    iter,
    alpha,
    psi,
    lo = NA_real_,
    hi = NA_real_,
    cached = FALSE
  ) {
    if (!verbose || iter %% trace_every != 0) {
      return(invisible(NULL))
    }

    msg <- sprintf(
      "%s iter=%03d alpha=%.8g psi=%.8g |psi-theta|=%.3g",
      stage,
      iter,
      alpha,
      psi,
      abs(psi - theta)
    )

    if (is.finite(lo) && is.finite(hi)) {
      msg <- paste0(
        msg,
        sprintf(" bracket=[%.8g, %.8g] width=%.3g", lo, hi, hi - lo)
      )
    }

    if (cached) {
      msg <- paste0(msg, " [cache hit]")
    }

    cat(msg, "\n")
    if (interactive()) {
      flush.console()
    }
    invisible(NULL)
  }

  # Initial evaluation
  rec0 <- eval_psi(alpha_init)
  record_eval(rec0)
  progress("init", 0L, rec0$alpha, rec0$psi, cached = rec0$cached)

  if (best_dist <= c_n) {
    ord <- order(alpha_hist[seq_len(n_eval)])
    return(list(
      alpha_hat = alpha_hist[best_idx],
      converged = TRUE,
      dist = best_dist,
      c_n = c_n,
      grid = data.frame(
        alpha = alpha_hist[seq_len(n_eval)][ord],
        psi = psi_hist[seq_len(n_eval)][ord],
        se = se_hist[seq_len(n_eval)][ord]
      ),
      eic = eic_hist[[best_idx]]
    ))
  }

  # This assumes psi(alpha) increases with alpha (Lemma 1 of Rytgaard and
  # van der Laan).
  if (rec0$psi < theta) {
    # Search upward to find an upper bracket
    lo <- rec0$alpha
    f_lo <- rec0$psi
    hi <- NA_real_
    f_hi <- NA_real_

    alpha <- rec0$alpha
    for (iter in seq_len(max_iter)) {
      alpha_new <- min(alpha_max, max(alpha_min, alpha * expand_up))
      if (alpha_new <= alpha * (1 + 1e-15)) {
        break
      }

      rec <- try(eval_psi(alpha_new), silent = TRUE)

      if (inherits(rec, "try-error")) {
        expand_up <- expand_up * 0.9
      } else {
        if (!isTRUE(rec$cached)) {
          record_eval(rec)
        }
        progress("expand\u2191", iter, rec$alpha, rec$psi, cached = rec$cached)

        if (rec$psi >= theta) {
          hi <- rec$alpha
          f_hi <- rec$psi
          break
        }

        lo <- rec$alpha
        f_lo <- rec$psi
        alpha <- rec$alpha
      }
    }
  } else {
    # Search downward to find a lower bracket
    hi <- rec0$alpha
    f_hi <- rec0$psi
    lo <- NA_real_
    f_lo <- NA_real_

    alpha <- rec0$alpha
    for (iter in seq_len(max_iter)) {
      alpha_new <- max(alpha_min, min(alpha_max, alpha * expand_down))
      if (alpha_new >= alpha * (1 - 1e-15)) {
        break
      }

      rec <- eval_psi(alpha_new)
      if (!isTRUE(rec$cached)) {
        record_eval(rec)
      }
      progress("expand\u2193", iter, rec$alpha, rec$psi, cached = rec$cached)

      if (rec$psi <= theta) {
        lo <- rec$alpha
        f_lo <- rec$psi
        break
      }

      hi <- rec$alpha
      f_hi <- rec$psi
      alpha <- rec$alpha
    }
  }

  bracketed <- is.finite(lo) &&
    is.finite(hi) &&
    lo < hi &&
    is.finite(f_lo) &&
    is.finite(f_hi) &&
    f_lo <= theta &&
    theta <= f_hi

  if (bracketed) {
    for (iter in seq_len(max_iter)) {
      if (abs(f_hi - f_lo) < .Machine$double.eps) {
        break
      }

      # Secant step, clipped to the bracket; fallback to midpoint if needed
      alpha_new <- lo + (theta - f_lo) * (hi - lo) / (f_hi - f_lo)
      if (!is.finite(alpha_new) || alpha_new <= lo || alpha_new >= hi) {
        alpha_new <- 0.5 * (lo + hi)
      }

      if (alpha_new <= alpha_min) {
        alpha_new <- alpha_min
      }
      if (alpha_new >= alpha_max) {
        alpha_new <- alpha_max
      }

      rec <- eval_psi(alpha_new)
      if (!isTRUE(rec$cached)) {
        record_eval(rec)
      }
      progress("search", iter, rec$alpha, rec$psi, lo, hi, cached = rec$cached)

      if (abs(rec$psi - theta) <= c_n) {
        break
      }

      if (rec$psi < theta) {
        lo <- rec$alpha
        f_lo <- rec$psi
      } else {
        hi <- rec$alpha
        f_hi <- rec$psi
      }

      if (abs(hi - lo) <= .Machine$double.eps * max(1, abs(lo), abs(hi))) break
    }
  } else if (verbose) {
    cat("No bracket found; returning the best point seen.\n")
    if (interactive()) flush.console()
  }

  # Best point seen
  alpha_best <- alpha_hist[best_idx]
  converged <- best_dist <= c_n

  ord <- order(alpha_hist[seq_len(n_eval)])

  list(
    alpha_hat = alpha_best,
    converged = converged,
    dist = best_dist,
    c_n = c_n,
    grid = data.frame(
      alpha = alpha_hist[seq_len(n_eval)][ord],
      psi = psi_hist[seq_len(n_eval)][ord],
      se = se_hist[seq_len(n_eval)][ord]
    ),
    eic = eic_hist[[best_idx]]
  )
}
