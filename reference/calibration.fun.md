# Calibrate a shape-parameter intervention to a target effect size

Finds the intensity-shape intervention `alpha` (via
[`estimate.alpha.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/estimate.alpha.fun.md)
and
[`tmle.alpha.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle.alpha.fun.md))
for the `z` process such that the resulting effect on `target` matches a
prespecified fraction `rho` of the effect under `alpha = 1`, or a fixed
`theta`/`delta`.

## Usage

``` r
calibration.fun(
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
)
```

## Arguments

- initial.fit:

  initial fit object, as returned by
  [`prepare.initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/prepare.initial.md).

- a:

  optional intervention on the baseline treatment `A0`, passed on to
  [`tmle.alpha.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle.alpha.fun.md).

- theta:

  optional numeric; target value for the estimand to calibrate to
  directly.

- rho:

  optional numeric; target fraction of the `alpha = 1` effect on
  `target` to calibrate to.

- delta:

  optional numeric; target absolute difference from the `alpha = 1`
  effect on `target` to calibrate to.

- rho.1a:

  optional numeric; as `rho`, but relative to the effect under `a` fixed
  to 1 rather than `alpha = 1`.

- browse:

  logical; if `TRUE`, drop into
  [`browser()`](https://rdrr.io/r/base/browser.html).

- verbose:

  logical; if `TRUE`, print progress.

- output.eic:

  logical; if `TRUE`, include influence-function values in the output.

- tau:

  follow-up horizon for the `target` process.

- tau.z:

  follow-up horizon for the `z` process (defaults to `tau`).

- use.cores:

  number of cores to use.

- target:

  character; name of the outcome process being calibrated to.

- z.name:

  character; name of the process on which the shape intervention `alpha`
  acts.

- ...:

  additional arguments passed on to
  [`tmle.alpha.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle.alpha.fun.md).

## Value

A list with the calibrated `alpha` (as found by
[`estimate.alpha.fun()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/estimate.alpha.fun.md))
and the corresponding TMLE fit(s).

## Examples

``` r
set.seed(1405)
baseline <- list(L0 = function(N) rnorm(N))
processes <- list(
  z = list(type = "one.jump", eta = 0.2, nu = 1),
  outcome1 = list(type = "terminal", eta = 0.3, nu = 1),
  censoring = list(type = "censoring", eta = 0.1, nu = 1)
)
effects <- list(
  c("L0", "z", 0.5),
  c("z", "outcome1", 0.7)
)
dt <- sim.generic(baseline = baseline, processes = processes, effects = effects, n = 30)

initial.fit <- prepare.initial(
  dt,
  tau = 1,
  fit.types = list(
    z = list(
      model = "Surv(tstart, tstop, delta == 2)~L0",
      fit = "cox",
      at.risk = function(dt) (dt[["z"]] == 0)
    ),
    outcome1 = list(model = "Surv(tstart, tstop, delta == 1)~L0+z", fit = "cox"),
    censoring = list(model = "Surv(tstart, tstop, delta == 0)~L0", fit = "cox")
  ),
  verbose = FALSE
)
#> Error in mclapply(state.chunks, function(state.chunk) {    chunk.out <- data.table(row_id = tmp.state[["row_id"]])    chunk.out[["state"]] <- 0    for (state.jj in state.chunk) {        if (verbose & (state.jj %in% print.state.jj)) {            print(paste0(state.jj, "/", max(depend.matrix[["state"]])))        }        depend.matrix.jj <- depend.matrix[state == state.jj]        which.jj <- rep(TRUE, nrow(tmp.state))        for (varname in setdiff(names(depend.matrix), "state")) {            tmp.state[[varname]] <- depend.matrix.jj[[varname]]            if (!any.hal) {                if (length(which.recurrent) > 0) {                  if (varname %in% which.recurrent) {                    for (count.jj in 1:max.count) {                      tmp.state[, `:=`((paste0(varname, ".",                         count.jj)), (tmp.state[[varname]] >=                         count.jj))]                    }                  }                }            }            which.jj <- which.jj & (tmp.long[[varname]] >= tmp.state[[varname]])        }        if (any.hal) {            X.hal.dynamic <- Matrix(model.matrix(formula(paste0("delta",                 "~-1+", paste(paste0("(", gsub("FALSE|TRUE",                   "", gsub(":", "):(", hal.vars.dynamic)), ")"),                   collapse = "+"))), data = tmp.state), sparse = TRUE)            X.hal.dynamic <- flip.interactions(X.hal.dynamic,                 hal.vars.dynamic)            if (length(hal.vars.A0.dynamic) > 0) {                X.hal.A0.dynamic <- Matrix(model.matrix(formula(paste0("delta",                   "~-1+", paste(paste0("(", gsub("FALSE|TRUE",                     "", gsub(":", "):(", hal.vars.A0.dynamic)),                     ")"), collapse = "+"))), data = tmp.state),                   sparse = TRUE)                X.hal.A0.dynamic <- flip.interactions(X.hal.A0.dynamic,                   hal.vars.A0.dynamic)            }            X.full <- cbind(X.hal.static, X.hal.dynamic, X.hal.A0.dynamic,                 X.hal.A0)        }        for (fit.type.jj in (1:length(fit.types))[-cens.process.id]) {            if (intervene.A0) {                colname <- paste0("P.a", aa, ".", names(fit.types)[fit.type.jj],                   ".", state.jj)            }            else {                colname <- paste0("P.", names(fit.types)[fit.type.jj],                   ".", state.jj)            }            if (fit.types[[fit.type.jj]]$fit == "hal") {                hal.index <- (1:length(which.hal))[names(which.hal) ==                   names(fit.types)[fit.type.jj]]                chunk.out[[colname]] <- as.numeric(exp(predict(fit.hals[[hal.index]][["hal.fit"]],                   X.full[, hal.vars.list[[hal.index]]], newoffset = 0,                   s = fit.hals[[hal.index]][["lambda.cv"]])) *                   tmp.state[["risk.time"]])            }            else {                lp <- predict(fit.cox.types[[fit.type.jj]]["fit.cox"][[1]],                   newdata = tmp.state, type = "lp")                chunk.out[[colname]] <- as.numeric(tmp.long[[paste0("dhazard.",                   names(fit.types)[fit.type.jj])]] * exp(lp))            }            if (fit.type.jj %in% at.risk.ids) {                chunk.out[[colname]] <- chunk.out[[colname]] *                   fit.types[[fit.type.jj]][["at.risk"]](depend.matrix.jj)            }        }        chunk.out[["state"]][which.jj == 1] <- state.jj    }    return(chunk.out)}, mc.cores = min(detectCores() - 5, use.cores.prediction, use.cores)): 'mc.cores' must be >= 1

# Calibrate alpha so the risk of outcome1 by tau = 1 equals 0.2.
cal <- calibration.fun(
  initial.fit = initial.fit,
  theta = 0.2,
  tau = 1,
  target = "outcome1",
  z.name = "z",
  use.cores = 1,
  verbose = FALSE
)
#> Error: object 'initial.fit' not found
cal$estimate
#> Error: object 'cal' not found
```
