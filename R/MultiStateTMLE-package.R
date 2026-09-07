#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import data.table
#' @importFrom Matrix Matrix
#' @importFrom nleqslv nleqslv
#' @importFrom parallel detectCores mclapply
#' @importFrom simevent simEventData
#' @importFrom stats as.formula binomial coef formula glm lm median
#'   model.matrix predict quantile sd setNames
#' @importFrom stringr str_split
#' @importFrom survival basehaz coxph coxph.control
#' @importFrom utils combn flush.console
#' @importFrom zoo na.locf
## usethis namespace: end
NULL

# Column names used via data.table's non-standard evaluation (e.g. `dt[, .(x)]`),
# plus a couple of arguments whose defaults reference an object expected to
# exist in the caller's environment. Not undefined globals at runtime; listed
# here only to silence `R CMD check`'s static analysis.
utils::globalVariables(c(
  "..stateful_derived",
  "..stateful_in_states",
  "A0.obs",
  "alpha.est.fun",
  "at.risk",
  "C",
  "C.1",
  "clever.weight",
  "clever.weight.alpha",
  "cum.hazard.z",
  "cum.hazard.z.1",
  "delta",
  "dt",
  "final.time",
  "grid.period",
  "grid.time",
  "Haz",
  "hazard",
  "id",
  "idN",
  "pi.A0.1",
  "Q",
  "risk.time",
  "row_id",
  "state",
  "state.row.index",
  "surv.0",
  "surv.0.1",
  "time",
  "time.obs",
  "tstart",
  "tstop",
  "x"
))
