# Targeted maximum likelihood estimation under a shape-parameter intervention

Estimates the effect of a stochastic intensity intervention (scaling the
shape of the `z` process by `alpha`, optionally combined with an
intervention `a` on baseline treatment) on `target`, via an iterative
TMLE update of the initial fit produced by
[`prepare.initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/prepare.initial.md).

## Usage

``` r
tmle.alpha.fun(
  target = "z",
  tau = 1.2,
  alpha = 1,
  z.name = "z",
  alpha.list = NULL,
  a = NULL,
  initial.fit = NULL,
  dt = NULL,
  years.lost = NULL,
  only.first = NULL,
  target.only.in.state = NULL,
  target.by.state = FALSE,
  conv.const = 1,
  one.step = FALSE,
  verbose = FALSE,
  max.iter = 10,
  min.iter = 0,
  use.cores = 50,
  truncate.weights = 0,
  output.convergence = FALSE,
  output.eic = FALSE,
  output.weights = NULL,
  output.a.weights = NULL,
  browse = FALSE,
  verbose.exponential = FALSE,
  ...
)
```

## Arguments

- target:

  character; name of the target process/outcome.

- tau:

  follow-up horizon (in time units) for `target`.

- alpha:

  numeric; multiplicative shape-parameter intervention applied to the
  `z` process (ignored if `alpha.list` is supplied).

- z.name:

  character; name of the process on which `alpha` acts.

- alpha.list:

  optional named list of intervention functions, one per process, of the
  form `function(time, covariates)`; supersedes `alpha`.

- a:

  optional intervention on the baseline treatment `A0` (must match the
  `a` used to build `initial.fit`, if `initial.fit` is supplied).

- initial.fit:

  initial fit object, as returned by
  [`prepare.initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/prepare.initial.md);
  computed from `dt`/`tau`/`a`/`...` if not supplied.

- dt:

  data.table of observed event history data, passed to
  [`prepare.initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/prepare.initial.md)
  when `initial.fit` is not supplied.

- years.lost:

  optional; if supplied, computes years-lost-type estimands with this
  block size.

- only.first:

  logical; if `TRUE`, restrict `target` to its first occurrence (treat
  it as a one-jump process).

- target.only.in.state:

  optional function of `states` identifying the subset of states in
  which `target` events are counted.

- target.by.state:

  logical; whether to compute clever covariates by state rather than
  pooled.

- conv.const:

  constant scaling the convergence criterion on the mean efficient
  influence curve.

- one.step:

  logical; if `TRUE`, return a one-step estimator instead of iterating
  to TMLE convergence.

- verbose:

  logical; if `TRUE`, print iteration progress.

- max.iter:

  maximum number of TMLE update iterations.

- min.iter:

  minimum number of TMLE update iterations before checking convergence.

- use.cores:

  number of cores to use for the per-id computations.

- truncate.weights:

  numeric; if `> 0`, truncate clever weights above this value.

- output.convergence:

  logical; if `TRUE`, include convergence diagnostics in the output.

- output.eic:

  logical; if `TRUE`, include the efficient influence curve in the
  output.

- output.weights:

  optional numeric vector of quantiles at which to summarize the
  clever/censoring weights.

- output.a.weights:

  optional numeric vector of quantiles at which to summarize the
  `a`-specific clever weights.

- browse:

  logical; if `TRUE`, drop into
  [`browser()`](https://rdrr.io/r/base/browser.html).

- verbose.exponential:

  logical; if `TRUE`, print diagnostics from the exponential TMLE update
  step.

- ...:

  additional arguments passed on to
  [`prepare.initial()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/prepare.initial.md)
  when `initial.fit` is not supplied.

## Value

A list with the TMLE `estimate`, standard error, and (depending on the
arguments above) the efficient influence curve, convergence diagnostics,
and weight summaries.

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
dt <- sim.generic(baseline = baseline, processes = processes, effects = effects, n = 100)

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

# Estimated risk of outcome1 by tau = 1 with no intervention (alpha = 1).
fit <- tmle.alpha.fun(
  initial.fit = initial.fit,
  target = "outcome1",
  tau = 1,
  alpha = 1,
  use.cores = 1,
  verbose = FALSE
)
#> Error: object 'initial.fit' not found
fit$estimate
#> Error: object 'fit' not found
```
