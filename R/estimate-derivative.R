#' Estimate the derivative of an estimate with respect to alpha
#'
#' Approximates the derivative with respect to `alpha` of the estimate
#' returned by `fun`, by a symmetric finite difference at `alpha`.
#'
#' @param fun A function of `alpha` (and `...`) returning a list with an
#'   element `estimate`, a named vector containing `tmle_est` or
#'   `one_step_est`, such as a call to [tmle_alpha()].
#' @param alpha Positive number. The point at which to estimate the
#'   derivative.
#' @param h Positive number. Step size of the finite difference. It is
#'   shrunk to `3 * alpha / 4` if it exceeds `alpha`, so that `fun` is only
#'   evaluated at positive values.
#' @param ... Further arguments passed on to `fun`.
#'
#' @return A number, the estimated derivative.
#'
#' @examples
#' # The derivative of alpha^2 at alpha = 3 is 6.
#' fun <- function(alpha) list(estimate = c(tmle_est = alpha^2))
#' estimate_derivative(fun, alpha = 3, h = 0.01)
#' @export
estimate_derivative <- function(fun, alpha, h, ...) {
  checkmate::assert_function(fun)
  checkmate::assert_number(alpha, lower = 0, finite = TRUE)
  checkmate::assert_number(h, lower = 0, finite = TRUE)

  if (alpha < h) {
    h <- alpha * 3 / 4
  }

  psi_plus <- fun(alpha = alpha + h, ...)
  psi_minus <- fun(alpha = alpha - h, ...)

  (.point_estimate(psi_plus) - .point_estimate(psi_minus)) / (2 * h)
}
