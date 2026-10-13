#' Prepare the initial fit (deprecated)
#'
#' `prepare.initial()` is deprecated as of MultiStateTMLE 0.1.0. Use
#' [fit_initial()] instead.
#'
#' @param dt Observed data, see `data` in [fit_initial()].
#' @param tau,a,verbose See [fit_initial()].
#' @param fit.types,max.count,derived.vars,depend.time,prune.states,use.exponential
#'   See the snake_case equivalents in [fit_initial()]. Use `at.risk` for the
#'   at-risk functions of `fit.types`.
#' @param fit.treatment See `treatment_model` in [fit_initial()].
#' @param cut.time,cut.one.way,cut.Tk,cut.Tk.values,two.way,reduce.NK,hal.sl,screen.two.way,lambda.cvs,event.dependent.cv,npenalize.vars,V,seed.hal,reduce.seed.dependence,penalize.time,verbose.hal
#'   HAL settings, see `hal_control` in [fit_initial()].
#' @param use.cores See `cores` in [fit_initial()].
#' @param verbose.exponential Merged into `verbose`.
#' @param use.cores.prediction,max.Tk,cv.glmnet,browse,browse.hal Ignored.
#' @param return.parameters.for.simulation No longer supported, since
#'   `sim.from.data()` has been removed. Use
#'   `simevent::sim_model_from_fits()` instead.
#'
#' @return See [fit_initial()].
#' @keywords internal
#' @export
prepare.initial <- function(
  dt,
  tau = 1.2,
  fit.types = list(
    Z = list(
      model = "Surv(tstart, tstop, delta == 2)~L0+L+A0",
      fit = "cox",
      at.risk = function(dt) (dt[["Z"]] == 0)
    ),
    L = list(
      model = "Surv(tstart, tstop, delta == 3)~L0+Z+A0",
      fit = "cox",
      at.risk = function(dt) (dt[["L"]] == 0)
    ),
    outcome1 = list(
      model = "Surv(tstart, tstop, delta == 1)~L0+Z+L+A0",
      fit = "cox"
    ),
    censoring = list(
      model = "Surv(tstart, tstop, delta == 0)~L0+Z+L+A0",
      fit = "cox"
    )
  ),
  use.exponential = FALSE,
  verbose.exponential = FALSE,
  browse = FALSE,
  max.count = 1,
  a = NULL,
  derived.vars = list(),
  depend.time = list(),
  fit.treatment = NULL,
  prune.states = FALSE,
  cut.time = 35,
  cut.one.way = 15,
  cut.Tk = 3,
  cut.Tk.values = NULL,
  max.Tk = max.count,
  two.way = NULL,
  reduce.NK = 0.025,
  hal.sl = list(c(
    cut.one.way = cut.one.way,
    cut.Tk = cut.Tk,
    cut.Tk.values = cut.Tk.values,
    two.way = two.way
  )),
  screen.two.way = FALSE,
  lambda.cvs = c(
    (9:1) / 10,
    (9:2) / 10^2,
    seq(1 / 10^2, 1 / 10^5, length = 100)
  ),
  event.dependent.cv = FALSE,
  npenalize.vars = NULL,
  V = 10,
  seed.hal = NULL,
  reduce.seed.dependence = FALSE,
  penalize.time = FALSE,
  use.cores = 1,
  use.cores.prediction = 1,
  verbose.hal = FALSE,
  browse.hal = FALSE,
  cv.glmnet = FALSE,
  verbose = FALSE,
  return.parameters.for.simulation = FALSE
) {
  lifecycle::deprecate_warn("0.1.0", "prepare.initial()", "fit_initial()")
  .prepare_initial(
    dt = dt,
    tau = tau,
    fit.types = fit.types,
    use.exponential = use.exponential,
    verbose.exponential = verbose.exponential,
    max.count = max.count,
    a = a,
    derived.vars = derived.vars,
    depend.time = depend.time,
    fit.treatment = fit.treatment,
    prune.states = prune.states,
    cut.time = cut.time,
    cut.one.way = cut.one.way,
    cut.Tk = cut.Tk,
    cut.Tk.values = cut.Tk.values,
    two.way = two.way,
    reduce.NK = reduce.NK,
    hal.sl = hal.sl,
    screen.two.way = screen.two.way,
    lambda.cvs = lambda.cvs,
    event.dependent.cv = event.dependent.cv,
    npenalize.vars = npenalize.vars,
    V = V,
    seed.hal = seed.hal,
    reduce.seed.dependence = reduce.seed.dependence,
    penalize.time = penalize.time,
    use.cores = use.cores,
    verbose.hal = verbose.hal,
    verbose = verbose,
    return.parameters.for.simulation = return.parameters.for.simulation
  )
}

# The body of prepare.initial(), without the deprecation warning, so that
# tmle.alpha.fun() can call it without warning twice. Arguments without a
# default here are always supplied by prepare.initial().
.prepare_initial <- function(
  dt,
  tau = 1.2,
  fit.types,
  use.exponential = FALSE,
  verbose.exponential = FALSE,
  max.count = 1,
  a = NULL,
  derived.vars = list(),
  depend.time = list(),
  fit.treatment = NULL,
  prune.states = FALSE,
  cut.time = 35,
  cut.one.way = 15,
  cut.Tk = 3,
  cut.Tk.values = NULL,
  two.way = NULL,
  reduce.NK = 0.025,
  hal.sl = NULL,
  screen.two.way = FALSE,
  lambda.cvs = .hal_control_defaults$lambda_cvs,
  event.dependent.cv = FALSE,
  npenalize.vars = NULL,
  V = 10,
  seed.hal = NULL,
  reduce.seed.dependence = FALSE,
  penalize.time = FALSE,
  use.cores = 1,
  verbose.hal = FALSE,
  verbose = FALSE,
  return.parameters.for.simulation = FALSE,
  ...
) {
  if (isTRUE(return.parameters.for.simulation)) {
    lifecycle::deprecate_stop(
      "0.1.0",
      "prepare.initial(return.parameters.for.simulation)",
      details = paste(
        "sim.from.data() has been removed;",
        "use simevent::sim_model_from_fits() instead."
      )
    )
  }
  fit_initial(
    data = dt,
    fit_types = fit.types,
    tau = tau,
    treatment_model = fit.treatment,
    a = a,
    max_count = max.count,
    derived_vars = derived.vars,
    depend_time = depend.time,
    prune_states = prune.states,
    use_exponential = use.exponential,
    hal_control = list(
      cut_time = cut.time,
      cut_one_way = cut.one.way,
      cut_tk = cut.Tk,
      cut_tk_values = cut.Tk.values,
      two_way = two.way,
      reduce_nk = reduce.NK,
      sl = hal.sl,
      screen_two_way = screen.two.way,
      lambda_cvs = lambda.cvs,
      event_dependent_cv = event.dependent.cv,
      npenalize_vars = npenalize.vars,
      v_folds = V,
      seed = seed.hal,
      reduce_seed_dependence = reduce.seed.dependence,
      penalize_time = penalize.time,
      verbose = verbose.hal
    ),
    cores = use.cores,
    verbose = verbose || verbose.exponential
  )
}

#' TMLE under a shape-parameter intervention (deprecated)
#'
#' `tmle.alpha.fun()` is deprecated as of MultiStateTMLE 0.1.0. Use
#' [tmle_alpha()] instead.
#'
#' @param initial.fit See `initial_fit` in [tmle_alpha()]. If `NULL`, it is
#'   computed by `prepare.initial(tau = tau, a = a, verbose = verbose, ...)`.
#' @param target,alpha,a,verbose See [tmle_alpha()].
#' @param z.name,alpha.list,years.lost,target.only.in.state,target.by.state,conv.const,one.step,max.iter,min.iter,truncate.weights,output.convergence,output.eic,output.weights,output.a.weights
#'   See the snake_case equivalents in [tmle_alpha()].
#' @param only.first See `only_first` in [tmle_alpha()]; `NULL` means
#'   `FALSE`.
#' @param use.cores See `cores` in [tmle_alpha()].
#' @param tau Only used when `initial.fit` is `NULL`; the horizon is that of
#'   `initial.fit`.
#' @param dt,browse,verbose.exponential Ignored.
#' @param ... Passed on to `prepare.initial()` when `initial.fit` is `NULL`.
#'
#' @return See [tmle_alpha()].
#' @keywords internal
#' @export
tmle.alpha.fun <- function(
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
  use.cores = 1,
  truncate.weights = 0,
  output.convergence = FALSE,
  output.eic = FALSE,
  output.weights = NULL,
  output.a.weights = NULL,
  browse = FALSE,
  verbose.exponential = FALSE,
  ...
) {
  lifecycle::deprecate_warn("0.1.0", "tmle.alpha.fun()", "tmle_alpha()")
  if (length(initial.fit) == 0) {
    initial.fit <- .prepare_initial(
      dt = dt,
      tau = tau,
      a = a,
      verbose = verbose,
      ...
    )
  }
  if (is.logical(output.weights)) {
    output.weights <- as.numeric(output.weights)
  }
  tmle_alpha(
    initial_fit = initial.fit,
    target = target,
    alpha = alpha,
    z_name = z.name,
    alpha_list = alpha.list,
    a = a,
    years_lost = years.lost,
    only_first = isTRUE(only.first),
    target_only_in_state = target.only.in.state,
    target_by_state = target.by.state,
    conv_const = conv.const,
    one_step = one.step,
    max_iter = max.iter,
    min_iter = min.iter,
    truncate_weights = truncate.weights,
    output_convergence = output.convergence,
    output_eic = output.eic,
    output_weights = output.weights,
    output_a_weights = output.a.weights,
    cores = use.cores,
    verbose = verbose
  )
}

#' Calibrate a shape-parameter intervention (deprecated)
#'
#' `calibration.fun()` is deprecated as of MultiStateTMLE 0.1.0. Use
#' [calibrate_alpha()] instead.
#'
#' @param initial.fit See `initial_fit` in [calibrate_alpha()].
#' @param target,a,theta,rho,delta,verbose See [calibrate_alpha()]. `theta`
#'   is ignored if `rho`, `delta` or `rho.1a` is given.
#' @param z.name,rho.1a,output.eic See the snake_case equivalents in
#'   [calibrate_alpha()].
#' @param use.cores See `cores` in [calibrate_alpha()].
#' @param tau,tau.z,browse Ignored; the horizon is that of `initial.fit`.
#' @param ... Passed on to [tmle_alpha()], with dotted argument names
#'   translated to their snake_case equivalents.
#'
#' @return See [calibrate_alpha()].
#' @keywords internal
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
  use.cores = 1,
  target = "outcome",
  z.name = "z",
  ...
) {
  lifecycle::deprecate_warn("0.1.0", "calibration.fun()", "calibrate_alpha()")
  if (length(c(rho, delta, rho.1a)) > 0) {
    theta <- NULL
  }
  do.call(
    calibrate_alpha,
    c(
      list(
        initial_fit = initial.fit,
        target = target,
        z_name = z.name,
        theta = theta,
        rho = rho,
        delta = delta,
        rho_1a = rho.1a,
        a = a,
        output_eic = output.eic,
        cores = use.cores,
        verbose = verbose
      ),
      .snake_case_args(list(...))
    )
  )
}

#' Calibration curve over a grid of interventions (deprecated)
#'
#' `calibration.curve.fun()` is deprecated as of MultiStateTMLE 0.1.0. Use
#' [calibration_curve()] instead.
#'
#' @param initial.fit See `initial_fit` in [calibration_curve()].
#' @param target,a,verbose See [calibration_curve()].
#' @param z.name,alpha.grid,min.iter,target.by.state,block.size.z,block.size.target,output.eic
#'   See the snake_case equivalents in [calibration_curve()].
#' @param use.cores See `cores` in [calibration_curve()].
#' @param tau,tau.z,browse Ignored; the horizon is that of `initial.fit`.
#' @param ... Passed on to [tmle_alpha()], with dotted argument names
#'   translated to their snake_case equivalents.
#'
#' @return See [calibration_curve()].
#' @keywords internal
#' @export
calibration.curve.fun <- function(
  initial.fit = NULL,
  a = NULL,
  alpha.grid = seq(0, 5, length = 10),
  browse = FALSE,
  verbose = TRUE,
  output.eic = FALSE,
  tau = 1.2,
  tau.z = tau,
  use.cores = 1,
  target = "outcome",
  z.name = "z",
  min.iter = 1,
  target.by.state = FALSE,
  block.size.z = NULL,
  block.size.target = NULL,
  ...
) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "calibration.curve.fun()",
    "calibration_curve()"
  )
  do.call(
    calibration_curve,
    c(
      list(
        initial_fit = initial.fit,
        target = target,
        z_name = z.name,
        alpha_grid = alpha.grid,
        a = a,
        min_iter = min.iter,
        target_by_state = target.by.state,
        block_size_z = block.size.z,
        block_size_target = block.size.target,
        output_eic = output.eic,
        cores = use.cores,
        verbose = verbose
      ),
      .snake_case_args(list(...))
    )
  )
}

#' Search for the alpha achieving a target value (deprecated)
#'
#' `estimate.alpha.fun()` is deprecated as of MultiStateTMLE 0.1.0. Use
#' [solve_alpha()] instead.
#'
#' @param theta,fun,c_n,alpha_init,expand_up,expand_down,max_iter,alpha_min,alpha_max,verbose,trace_every
#'   See [solve_alpha()].
#' @param use.cores Ignored.
#'
#' @return See [solve_alpha()].
#' @keywords internal
#' @export
estimate.alpha.fun <- function(
  theta,
  fun,
  c_n,
  alpha_init = 1,
  expand_up = 2,
  expand_down = 0.25,
  max_iter = 100,
  alpha_min = 1e-3,
  alpha_max = 100,
  use.cores = 1,
  verbose = FALSE,
  trace_every = 1L
) {
  lifecycle::deprecate_warn("0.1.0", "estimate.alpha.fun()", "solve_alpha()")
  solve_alpha(
    fun = fun,
    theta = theta,
    c_n = c_n,
    alpha_init = alpha_init,
    expand_up = expand_up,
    expand_down = expand_down,
    max_iter = max_iter,
    alpha_min = alpha_min,
    alpha_max = alpha_max,
    verbose = verbose,
    trace_every = trace_every
  )
}

#' Estimate the derivative with respect to alpha (deprecated)
#'
#' `estimate.derivative()` is deprecated as of MultiStateTMLE 0.1.0. Use
#' [estimate_derivative()] instead.
#'
#' @param alpha_hat See `alpha` in [estimate_derivative()].
#' @param fun,h See [estimate_derivative()].
#' @param parameter Passed on to `fun` as `parameter`.
#' @param ... Passed on to `fun`.
#'
#' @return See [estimate_derivative()].
#' @keywords internal
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
  lifecycle::deprecate_warn(
    "0.1.0",
    "estimate.derivative()",
    "estimate_derivative()"
  )
  estimate_derivative(fun, alpha = alpha_hat, h = h, parameter = parameter, ...)
}

#' Total, mediated and calibrated contrasts (deprecated)
#'
#' `make.calibrated.contrasts()` is deprecated as of MultiStateTMLE 0.1.0.
#' Use [calibrated_contrasts()] instead.
#'
#' @param target See [calibrated_contrasts()].
#' @param treatment.fit,placebo.fit,calibrated.fit,conf.level See the
#'   snake_case equivalents in [calibrated_contrasts()].
#'
#' @return See [calibrated_contrasts()].
#' @keywords internal
#' @export
make.calibrated.contrasts <- function(
  target = NULL,
  treatment.fit,
  placebo.fit,
  calibrated.fit,
  conf.level = 0.95
) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "make.calibrated.contrasts()",
    "calibrated_contrasts()"
  )
  calibrated_contrasts(
    treatment_fit = treatment.fit,
    placebo_fit = placebo.fit,
    calibrated_fit = calibrated.fit,
    target = target,
    conf_level = conf.level
  )
}

#' Compute Q and clever covariates per id (deprecated)
#'
#' `compute.Q.clever.per.id()` is deprecated as of MultiStateTMLE 0.1.0 and
#' is now internal; it is called by [tmle_alpha()].
#'
#' @param dt_id A `data.table` with the time-ordered rows of a single id,
#'   with hazard columns named `<P.prefix><name>.<s>` for each process
#'   `<name>` and state index `s`.
#' @param states A `data.table` enumerating the product state space, with a
#'   column per stateful process and the state index column.
#' @param process.types Named list of process types (`"terminal"`,
#'   `"recurrent"`, `"one.jump"`).
#' @param P.prefix Prefix of the hazard columns.
#' @param state.idx.col Name of the state index column.
#' @param parameter Name of the target process; `"target"` picks the first
#'   terminal process.
#' @param compute.clever Logical. Whether to compute clever covariates.
#' @param get.years.lost,years.lost.block.size Whether to compute the
#'   expected time lost, and the block size of its approximation.
#' @param clever.by.state Logical. Whether to compute state-specific clever
#'   covariates.
#' @param process.deltas,browse,browse2 Ignored.
#'
#' @return `dt_id` with `Q`, `Q.<s>` and `clever.Q.<name>` columns added.
#' @keywords internal
#' @export
compute.Q.clever.per.id <- function(
  dt_id,
  states,
  process.types = NULL,
  P.prefix = "P.",
  state.idx.col = "state",
  parameter = "target",
  process.deltas = NULL,
  compute.clever = TRUE,
  browse = FALSE,
  browse2 = FALSE,
  get.years.lost = FALSE,
  years.lost.block.size = 10,
  clever.by.state = FALSE
) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "compute.Q.clever.per.id()",
    details = "It is now internal to tmle_alpha()."
  )
  .compute_q_clever(
    dt_id = dt_id,
    states = states,
    process_types = process.types,
    p_prefix = P.prefix,
    state_idx_col = state.idx.col,
    parameter = parameter,
    compute_clever = compute.clever,
    get_years_lost = get.years.lost,
    years_lost_block_size = years.lost.block.size,
    clever_by_state = clever.by.state
  )
}

# Translate the dotted argument names that the deprecated functions passed
# on to tmle.alpha.fun() through `...` to those of tmle_alpha(), dropping the
# ones it no longer has.
.snake_case_args <- function(args) {
  if (length(args) == 0) {
    return(args)
  }
  nms <- names(args)
  nms[nms == "use.cores"] <- "cores"
  nms <- gsub(".", "_", nms, fixed = TRUE)
  names(args) <- nms
  args[!(nms %in% c("tau", "dt", "browse", "verbose_exponential"))]
}
