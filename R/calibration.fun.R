### calibration.fun.R ---
#----------------------------------------------------------------------
## Author: Helene
## Created: May 13 2026 (19:31)
## Version:
## Last-Updated: Aug 28 2026 (08:59)
##           By: Helene
##     Update #: 140
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:

#' Calibrate a shape-parameter intervention to a target effect size
#'
#' Finds the intensity-shape intervention `alpha` (via [estimate.alpha.fun()]
#' and [tmle.alpha.fun()]) for the `z` process such that the resulting effect
#' on `target` matches a prespecified fraction `rho` of the effect under
#' `alpha = 1`, or a fixed `theta`/`delta`.
#'
#' @param initial.fit initial fit object, as returned by [prepare.initial()].
#' @param a optional intervention on the baseline treatment `A0`, passed on to
#'   [tmle.alpha.fun()].
#' @param theta optional numeric; target value for the estimand to calibrate
#'   to directly.
#' @param rho optional numeric; target fraction of the `alpha = 1` effect on
#'   `target` to calibrate to.
#' @param delta optional numeric; target absolute difference from the
#'   `alpha = 1` effect on `target` to calibrate to.
#' @param rho.1a optional numeric; as `rho`, but relative to the effect under
#'   `a` fixed to 1 rather than `alpha = 1`.
#' @param browse logical; if `TRUE`, drop into `browser()`.
#' @param verbose logical; if `TRUE`, print progress.
#' @param output.eic logical; if `TRUE`, include influence-function values in
#'   the output.
#' @param tau follow-up horizon for the `target` process.
#' @param tau.z follow-up horizon for the `z` process (defaults to `tau`).
#' @param use.cores number of cores to use.
#' @param target character; name of the outcome process being calibrated to.
#' @param z.name character; name of the process on which the shape
#'   intervention `alpha` acts.
#' @param ... additional arguments passed on to [tmle.alpha.fun()].
#' @return A list with the calibrated `alpha` (as found by
#'   [estimate.alpha.fun()]) and the corresponding TMLE fit(s).
#' @examples
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
#' dt <- sim.generic(baseline = baseline, processes = processes, effects = effects, n = 30)
#'
#' initial.fit <- prepare.initial(
#'   dt,
#'   tau = 1,
#'   fit.types = list(
#'     z = list(
#'       model = "Surv(tstart, tstop, delta == 2)~L0",
#'       fit = "cox",
#'       at.risk = function(dt) (dt[["z"]] == 0)
#'     ),
#'     outcome1 = list(model = "Surv(tstart, tstop, delta == 1)~L0+z", fit = "cox"),
#'     censoring = list(model = "Surv(tstart, tstop, delta == 0)~L0", fit = "cox")
#'   ),
#'   verbose = FALSE
#' )
#'
#' # Calibrate alpha so the risk of outcome1 by tau = 1 equals 0.2.
#' cal <- calibration.fun(
#'   initial.fit = initial.fit,
#'   theta = 0.2,
#'   tau = 1,
#'   target = "outcome1",
#'   z.name = "z",
#'   use.cores = 1,
#'   verbose = FALSE
#' )
#' cal$estimate
#' @export
calibration.fun <- function(
  initial.fit = NULL,
  a = NULL,
  theta = 0.5,
  rho = NULL,
  delta = NULL,
  rho.1a = NULL,
  browse = FALSE,
  verbose = TRUE,
  output.eic = FALSE,
  tau = 1.2,
  tau.z = tau,
  use.cores = 50,
  target = "outcome",
  z.name = "z",
  ...
) {
  a.fixed <- a
  use.cores.fixed <- use.cores
  verbose.fixed <- verbose
  tau.fixed <- tau
  z.name.fixed <- z.name

  fit.local <- data.table::copy(initial.fit)

  n <- length(unique(initial.fit[["tmp.long"]][["id"]]))

  tmle.alpha.fun.fixed <- function(
    tau = tau.z,
    alpha = 1,
    parameter = z.name,
    output.eic = TRUE,
    a = a.fixed
  ) {
    tmle.alpha.fun(
      initial.fit = data.table::copy(fit.local),
      tau = tau,
      alpha = alpha,
      a = a,
      target = parameter,
      z.name = z.name.fixed,
      use.cores = use.cores.fixed,
      verbose = verbose.fixed,
      output.eic = output.eic,
      ...
    )
  }

  if (length(rho) > 0) {
    est.aux.1 <- tmle.alpha.fun.fixed(alpha = 1, a = a, output.eic = TRUE)
    theta <- rho *
      est.aux.1[["estimate"]][grep(
        "tmle.est|one.step.est",
        names(est.aux.1[["estimate"]]),
        value = TRUE
      )]
    theta.eic <- rho * est.aux.1$eic
  } else if (length(rho.1a) > 0) {
    est.aux.1 <- tmle.alpha.fun.fixed(alpha = 1, a = 1 - a, output.eic = TRUE)
    theta <- rho.1a *
      est.aux.1[["estimate"]][grep(
        "tmle.est|one.step.est",
        names(est.aux.1[["estimate"]]),
        value = TRUE
      )]
    theta.eic <- rho.1a * est.aux.1$eic
  } else if (length(delta) > 0) {
    est.aux.1 <- tmle.alpha.fun.fixed(alpha = 1, a = a, output.eic = TRUE)
    theta <- delta +
      est.aux.1[["estimate"]][grep(
        "tmle.est|one.step.est",
        names(est.aux.1[["estimate"]]),
        value = TRUE
      )]
    theta.eic <- est.aux.1$eic
  } else {
    theta.eic <- 0
  }

  theta.se <- sqrt(mean(theta.eic^2) / n)

  if (theta < 0) {
    stop(paste0("inadmissible target: theta = ", theta))
  }

  if (browse) {
    browser()
  }

  if (theta > 0) {
    est.alpha <- estimate.alpha.fun(
      fun = tmle.alpha.fun.fixed,
      c_n = n^{
        -1 / 2
      } /
        log(n),
      verbose = verbose,
      theta = theta
    )

    if (!isTRUE(est.alpha$converged)) {
      stop(paste0(
        "calibration failed to converge: no alpha achieves theta = ",
        theta,
        " (closest achieved: ",
        est.alpha$grid$psi[which.min(abs(est.alpha$grid$psi - theta))],
        ", dist = ",
        est.alpha$dist,
        "); theta may exceed L(P), the achievable ceiling of Psi_z^alpha(P)"
      ))
    }

    est.deriv.auxiliary <- estimate.derivative(
      est.alpha$alpha.hat,
      parameter = z.name,
      fun = tmle.alpha.fun.fixed,
      h = 0.3 *
        n^{
          -1 / 6
        } *
        sqrt(mean(est.alpha$eic^2))
    )

    alpha.eic <- (1 / est.deriv.auxiliary) * (theta.eic - est.alpha$eic)
  } else {
    est.alpha <- list(
      alpha.hat = 0,
      eic = 0,
      converged = TRUE,
      dist = 0,
      c.n = 0,
      grid = data.frame(alpha = 0, psi = 0)
    )
    alpha.eic <- 0
    est.deriv.auxiliary <- Inf
  }

  alpha.se <- sqrt(mean(alpha.eic^2) / n)

  target.est <- tmle.alpha.fun.fixed(
    tau = tau,
    parameter = target,
    alpha = est.alpha$alpha.hat,
    output.eic = TRUE
  )

  target.se.crude <- target.est[["estimate"]]["se"]

  if (theta > 0) {
    est.deriv.target <- estimate.derivative(
      est.alpha$alpha.hat,
      parameter = target,
      tau = tau,
      fun = tmle.alpha.fun.fixed,
      h = 0.3 *
        n^{
          -1 / 6
        } *
        sqrt(mean(est.alpha$eic^2))
    )

    target.se <- sqrt(
      mean((target.est[["eic"]] + est.deriv.target * alpha.eic)^2) / n
    )
  } else {
    est.deriv.target <- Inf
    target.se <- sqrt(mean((target.est[["eic"]])^2) / n)
  }

  out.estimate <- c(
    alpha.est = est.alpha$alpha.hat,
    alpha.se = alpha.se,
    target.est = as.numeric(target.est[["estimate"]][grep(
      "tmle.est|one.step.est",
      names(target.est[["estimate"]]),
      value = TRUE
    )]), #["tmle.est"]
    target.se = target.se,
    target.se.crude = as.numeric(target.se.crude),
    theta = theta,
    theta.se = theta.se
  )

  out.checks <- c(
    est.deriv.auxiliary = est.deriv.auxiliary,
    est.alpha.converged = est.alpha$converged,
    est.alpha.dist = est.alpha$dist,
    est.alpha.cn = est.alpha$c.n,
    est.deriv.target = est.deriv.target
  )

  if (output.eic) {
    return(list(
      estimate = out.estimate,
      checks = out.checks,
      alpha.eic = alpha.eic,
      target.eic = target.est[["eic"]],
      eic = target.est[["eic"]] + est.deriv.target * alpha.eic
    ))
  } else {
    return(list(estimate = out.estimate, checks = out.checks))
  }
}

######################################################################
### calibration.fun.R ends here
