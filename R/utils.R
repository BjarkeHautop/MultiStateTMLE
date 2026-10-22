# Defaults for the `control` argument of fit_initial().
.fit_control_defaults <- list(
  prune_states = FALSE,
  use_exponential = FALSE
)

# Defaults for the `control` argument of tmle_alpha().
.tmle_control_defaults <- list(
  one_step = FALSE,
  block_size = 10,
  target_by_state = FALSE,
  conv_const = 1,
  max_iter = 10,
  min_iter = 0,
  truncate_weights = Inf
)

# Validate the `control` argument of tmle_alpha() and fill in its defaults.
.tmle_control <- function(control) {
  checkmate::assert_list(control, names = "unique")
  checkmate::assert_subset(names(control), names(.tmle_control_defaults))
  control <- utils::modifyList(.tmle_control_defaults, control)
  checkmate::assert_flag(control$one_step, .var.name = "control$one_step")
  checkmate::assert_count(
    control$block_size,
    positive = TRUE,
    .var.name = "control$block_size"
  )
  checkmate::assert_flag(
    control$target_by_state,
    .var.name = "control$target_by_state"
  )
  checkmate::assert_number(
    control$conv_const,
    lower = .Machine$double.eps,
    finite = TRUE,
    .var.name = "control$conv_const"
  )
  checkmate::assert_count(
    control$max_iter,
    positive = TRUE,
    .var.name = "control$max_iter"
  )
  checkmate::assert_count(control$min_iter, .var.name = "control$min_iter")
  checkmate::assert_number(
    control$truncate_weights,
    lower = .Machine$double.eps,
    .var.name = "control$truncate_weights"
  )
  control
}

# Defaults for the `hal_control` argument of fit_initial().
.hal_control_defaults <- list(
  cut_time = 35,
  cut_one_way = 15,
  cut_tk = 3,
  cut_tk_values = NULL,
  two_way = NULL,
  reduce_nk = 0.025,
  sl = NULL,
  screen_two_way = FALSE,
  lambda_cvs = c(
    (9:1) / 10,
    (9:2) / 10^2,
    seq(1 / 10^2, 1 / 10^5, length = 100)
  ),
  event_dependent_cv = FALSE,
  npenalize_vars = NULL,
  v_folds = 10,
  seed = NULL,
  reduce_seed_dependence = FALSE,
  penalize_time = FALSE,
  verbose = FALSE
)

# Bring one entry of fit_initial()'s `fit_types`,
# list(delta, model, fit, at_risk), to the internal form
# list(delta, model, fit, at.risk), where `model` is the string
# "Surv(tstart,tstop,delta==<delta>)~<covariates>" without whitespace, as
# the string manipulations on it in fit_initial() expect.
.normalize_fit_type <- function(spec, name) {
  var_name <- paste0("fit_types[[\"", name, "\"]]")
  checkmate::assert_list(spec, names = "unique", .var.name = var_name)
  checkmate::assert_subset(
    names(spec),
    c("delta", "model", "fit", "at_risk"),
    .var.name = paste0("names(", var_name, ")")
  )
  checkmate::assert_int(
    spec$delta,
    lower = 0,
    .var.name = paste0(var_name, "$delta")
  )
  checkmate::assert_formula(spec$model, .var.name = paste0(var_name, "$model"))
  if (length(spec$model) != 2) {
    stop(
      var_name,
      "$model must be a one-sided formula of the covariates, e.g. ~ L0; ",
      "the response is built from `delta`."
    )
  }
  covariates <- .deparse_one_line(spec$model[[2]])
  spec$model <- gsub(
    "\\s+",
    "",
    paste0("Surv(tstart,tstop,delta==", spec$delta, ")~", covariates)
  )
  if (is.null(spec$fit)) {
    spec$fit <- "cox"
  }
  checkmate::assert_choice(
    spec$fit,
    c("cox", "hal"),
    .var.name = paste0(var_name, "$fit")
  )
  if (!is.null(spec$at_risk)) {
    checkmate::assert_function(
      spec$at_risk,
      .var.name = paste0(var_name, "$at_risk")
    )
    spec$at.risk <- spec$at_risk
    spec$at_risk <- NULL
  }
  spec
}

# Number of cores to hand to parallel::mclapply(), which does not support
# more than one on Windows.
.n_cores <- function(cores) {
  if (.Platform$OS.type == "windows") {
    return(1L)
  }
  as.integer(max(1L, min(cores, parallel::detectCores(), na.rm = TRUE)))
}

# A language object, such as a formula, as a single string.
.deparse_one_line <- function(x) {
  paste(deparse(x, width.cutoff = 500L), collapse = "")
}

# The estimate table of the result objects: one row per target with the
# estimate, its standard error and a 95% Wald confidence interval.
.wald_table <- function(target, estimate, se) {
  z <- stats::qnorm(0.975)
  data.table(
    target = target,
    estimate = estimate,
    se = se,
    lower = estimate - z * se,
    upper = estimate + z * se
  )
}

# The point estimate of a tmle_alpha() fit.
.point_estimate <- function(fit) {
  as.numeric(fit[["estimate"]][["estimate"]][1])
}

# The efficient influence curve of a tmle_alpha() fit.
.eic <- function(fit) {
  fit[["eic"]][[1]]
}
