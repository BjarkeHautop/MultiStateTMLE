### estimate.derivative.R ---
#----------------------------------------------------------------------
## Author: Helene
## Created: May 13 2026 (19:43)
## Version:
## Last-Updated: Jun 25 2026 (15:51)
##           By: Helene
##     Update #: 22
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:

#' Estimate the numerical derivative of an alpha-indexed estimand
#'
#' Approximates the derivative with respect to `alpha` of the TMLE/one-step
#' estimate returned by `fun`, using a symmetric finite-difference at
#' `alpha_hat`.
#'
#' @param alpha_hat numeric; the point at which to estimate the derivative.
#' @param fun function taking `alpha`, `parameter`, and `...`, returning a
#'   list with an `estimate` element containing `tmle.est`/`one.step.est`.
#' @param parameter character; the target process/outcome passed to `fun`.
#' @param h numeric step size for the finite difference; shrunk if it exceeds
#'   `alpha_hat`.
#' @param ... additional arguments passed on to `fun`.
#' @return numeric; the estimated derivative.
#' @export
estimate.derivative <- function(
  alpha_hat,
  fun = alpha.est.fun,
  parameter = "z",
  h = length(unique(dt[, id]))^{
    -1 / 6
  },
  ...
) {
  if (alpha_hat < h) {
    h <- alpha_hat * 3 / 4
  }

  psi_plus <- fun(alpha = alpha_hat + h, parameter = parameter, ...)

  psi_minus <- fun(alpha = alpha_hat - h, parameter = parameter, ...)

  return(
    (psi_plus[["estimate"]][grep(
      "tmle.est|one.step.est",
      names(psi_plus[["estimate"]]),
      value = TRUE
    )] -
      psi_minus[["estimate"]][grep(
        "tmle.est|one.step.est",
        names(psi_minus[["estimate"]]),
        value = TRUE
      )]) /
      (2 * h)
  )
}

######################################################################
### estimate.derivative.R ends here
