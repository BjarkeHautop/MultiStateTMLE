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
  covariates <- paste(
    deparse(spec$model[[2]], width.cutoff = 500L),
    collapse = ""
  )
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

# Name the per-target values `x` (named by target, or "" for a single target)
# as `prefix`, suffixed by `_<target>` when there are several targets.
.target_names <- function(x, prefix) {
  names(x) <- if (identical(names(x), "")) {
    prefix
  } else {
    paste0(prefix, "_", names(x))
  }
  x
}

# The point estimate of a tmle_alpha() fit (TMLE or one-step).
.point_estimate <- function(fit) {
  est <- fit[["estimate"]]
  nm <- intersect(c("tmle_est", "one_step_est"), names(est))
  as.numeric(est[[nm[1]]])
}
