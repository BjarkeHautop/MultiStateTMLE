#' Simulated illness-death trial data
#'
#' Simulated event history data for a trial of 200 individuals randomised
#' to treatment (`A0 = 1`) or placebo (`A0 = 0`), at risk of `illness`,
#' which can occur at most once, and `death`. Both have exponential hazards.
#' Treatment lowers the hazard of the illness, which raises the hazard of
#' death, and also lowers the hazard of death directly, so its effect on
#' death is partly through the illness. The illness hazard also depends on a
#' standard normal baseline covariate `L0`, and follow-up is subject to
#' independent censoring. Simulated with [simevent::simulate_events()], see
#' `data-raw/illness_death.R` in the package source.
#'
#' @format A data frame with one row per event and columns:
#' \describe{
#'   \item{id}{Individual identifier.}
#'   \item{time}{Event time.}
#'   \item{event}{The event: `"illness"`, `"death"` or
#'     `"censoring"`.}
#'   \item{L0}{Baseline covariate.}
#'   \item{A0}{Randomised baseline treatment, `1` for treatment and `0` for
#'     placebo.}
#' }
"illness_death"
