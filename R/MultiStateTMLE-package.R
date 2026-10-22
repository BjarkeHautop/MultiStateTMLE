#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import data.table
#' @importFrom Matrix Matrix
#' @importFrom nleqslv nleqslv
#' @importFrom parallel detectCores mclapply
#' @importFrom stats as.formula binomial coef confint formula glm lm median
#'   model.matrix predict quantile sd setNames
#' @importFrom survival basehaz coxph coxph.control
#' @importFrom utils combn flush.console
## usethis namespace: end
NULL

# Column names used via data.table's non-standard evaluation (e.g. `dt[, .(x)]`),
# and fit.hal(), which HAL fits expect on the search path. Listed here only to
# silence `R CMD check`'s static analysis.
utils::globalVariables(c(
  "A0.obs",
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
  "fit.hal",
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
  "target.only.in.state",
  "time",
  "time.obs",
  "tstart",
  "tstop",
  "x"
))
