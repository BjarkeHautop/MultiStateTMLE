### sim.generic.R ---
#----------------------------------------------------------------------
## Author: Helene
## Created: Aug 28 2026 (09:46)
## Version:
## Last-Updated: Aug 30 2026 (19:50)
##           By: Helene
##     Update #: 290
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:

#' Simulate event history data from a generic process specification
#'
#' Simulates multistate event history data from a set of user-specified
#' baseline covariates, event processes (with Weibull intensities and Cox-type
#' effects), and their effects on one another, via [simevent::simEventData()].
#'
#' @param baseline named list of baseline covariate generators.
#' @param processes named list of process specifications; each entry has at
#'   least a `type` (one of `"censoring"`, `"terminal"`, `"one.jump"`, or
#'   recurrent), and Weibull intensity parameters `eta`/`nu`.
#' @param effects list of `c(from, to, coefficient)` triples specifying Cox
#'   effects between baseline/process variables and process intensities.
#' @param sim.object optional list with `baseline`/`processes`/`effects`
#'   entries, used as a fallback when none of `baseline`, `processes`, or
#'   `effects` are supplied directly.
#' @param cens numeric at-risk indicator scaling for the censoring process.
#' @param alpha.intervention named list of multiplicative interventions on
#'   process intensities (`eta`), keyed by process name.
#' @param baseline.intervention named list of interventions that fix a
#'   baseline covariate to a constant value.
#' @param n number of individuals to simulate.
#' @param browse logical; if `TRUE`, drop into `browser()` before simulating.
#' @return A `data.table` of simulated event history data with columns `id`,
#'   `time`, `delta`, baseline covariates, and one column per process.
#' @examples
#' # An illness-death process: z is a one-jump "illness" event that raises
#' # the hazard of the terminal process outcome1, with independent censoring.
#' set.seed(1405)
#' baseline <- list(L0 = function(N) rnorm(N))
#' processes <- list(
#'   z = list(type = "one.jump", eta = 0.2, nu = 1),
#'   outcome1 = list(type = "terminal", eta = 0.3, nu = 1),
#'   censoring = list(type = "censoring", eta = 0.1, nu = 1)
#' )
#' effects <- list(
#'   c("L0", "z", 0.5),
#'   c("z", "outcome1", 0.7)
#' )
#'
#' dt <- sim.generic(baseline = baseline, processes = processes, effects = effects, n = 100)
#' head(dt)
#' @export
sim.generic <- function(
  baseline = list(),
  processes = list(),
  effects = list(),
  sim.object = list(),
  cens = 1,
  alpha.intervention = list(),
  baseline.intervention = list(),
  n = 500,
  browse = FALSE
) {
  if (missing(baseline) && missing(processes) && missing(effects)) {
    baseline <- sim.object$baseline
    processes <- sim.object$processes
    effects <- sim.object$effects
  }

  baseline.vars <- names(baseline)

  add_cov <- copy(baseline)

  process.names <- names(processes)

  which.cens <- process.names[sapply(processes, function(process) {
    process[["type"]] == "censoring"
  })]
  which.terminal <- process.names[sapply(processes, function(process) {
    process[["type"]] == "terminal"
  })]
  which.one.jump <- setdiff(
    process.names[sapply(processes, function(process) {
      process[["type"]] == "one.jump"
    })],
    c(which.terminal, which.cens)
  )

  process.order <- c(
    which.cens,
    which.terminal,
    setdiff(process.names, c(which.cens, which.terminal))
  )

  eta <- sapply(processes[process.order], function(process) process[["eta"]])
  nu <- sapply(processes[process.order], function(process) process[["nu"]])

  if (length(baseline.intervention) > 0) {
    for (bname in names(baseline.intervention)) {
      add_cov[[bname]] <- function(N) rep(baseline.intervention[[bname]], N)
    }
  }

  if (length(alpha.intervention) > 0) {
    for (alphaname in names(alpha.intervention)) {
      eta[names(processes[process.order]) == alphaname] <-
        alpha.intervention[[alphaname]] *
        eta[names(processes[process.order]) == alphaname]
    }
  }

  at_risk <- function(events) {
    out <- numeric(length(process.order))

    names(out) <- process.order

    ## censoring
    out[process.order %in% which.cens] <- cens

    ## terminal events
    out[process.order %in% which.terminal] <- 1

    ## one jump
    for (one.jump in which.one.jump) {
      idx <- which(process.order == one.jump)
      out[idx] <- as.numeric(events[idx] == 0)
    }

    ## recurrent
    out[setdiff(
      process.order,
      c(which.cens, which.terminal, which.one.jump)
    )] <- 1

    return(out)
  }

  if (!("A0" %in% baseline.vars)) {
    add_A0 <- 1
  } else {
    add_A0 <- 0
  }

  if (!("L0" %in% baseline.vars)) {
    add_L0 <- 1
  } else {
    add_L0 <- 0
  }

  other.baseline.vars <- setdiff(baseline.vars, c("L0", "A0"))

  beta <- matrix(
    0,
    nrow = length(process.order) + length(other.baseline.vars) + 2,
    ncol = length(process.order)
  )

  # simEventData() always renames beta's rows to L0, A0, ... positionally
  # (to match its internal simmatrix), so this order must be fixed
  # regardless of whether L0/A0 were user-supplied or auto-added.
  rownames(beta) <- c("L0", "A0", other.baseline.vars, process.order)
  colnames(beta) <- process.order

  for (effect in effects) {
    beta[effect[1], effect[2]] <- as.numeric(effect[3])
  }

  override_beta <- NULL

  if (browse) {
    browser()
  }

  term.processes <- c(which.cens, which.terminal)
  term.deltas <- match(term.processes, process.order) - 1L

  non.term.processes <- setdiff(process.order, term.processes)
  non.term.deltas <- match(non.term.processes, process.order) - 1L

  data <- simEventData(
    N = n,
    beta = beta,
    eta = eta,
    nu = nu,
    max_cens = Inf,
    max_events = 50,
    at_risk = at_risk,
    lower = 1e-25,
    upper = 1e8,
    term_deltas = term.deltas,
    gen_L0 = add_cov[["L0"]],
    gen_A0 = {
      if ("A0" %in% names(add_cov)) function(N, L0) add_cov[["A0"]](N) else NULL
    },
    add_cov = add_cov[!(names(add_cov) %in% c("A0", "L0"))],
    override_beta = override_beta
  )

  if (add_L0) {
    data[["L0"]] <- NULL
  }

  if (add_A0) {
    data[["A0"]] <- NULL
  }

  for (jj in term.deltas) {
    data[[paste0("N", jj)]] <- NULL
  }

  if (length(non.term.processes) > 0) {
    setnames(data, paste0("N", non.term.deltas), non.term.processes)
  }

  setnames(data, c("Delta", "Time", "ID"), c("delta", "time", "id"))

  return(data)
}

######################################################################
### sim.generic.R ends here
