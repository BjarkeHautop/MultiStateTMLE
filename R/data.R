#' Simulated illness-death data
#'
#' Simulated event history data for 150 individuals at risk of an illness
#' `z`, which can occur at most once, and of death `outcome1`. Both have
#' exponential hazards that depend on a standard normal baseline covariate
#' `L0`, the illness raises the hazard of death, and follow-up is subject to
#' independent censoring. Simulated with [simevent::sim_events()], see
#' `data-raw/illness_death.R` in the package source.
#'
#' @format A data frame with one row per event and columns:
#' \describe{
#'   \item{id}{Individual identifier.}
#'   \item{time}{Event time.}
#'   \item{delta}{Event type: `0` censoring, `1` death (`outcome1`), `2`
#'     illness (`z`).}
#'   \item{L0}{Baseline covariate.}
#' }
"illness_death"
