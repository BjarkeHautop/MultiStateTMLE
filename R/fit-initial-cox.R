# Cox hazard models for fit_initial().

# Fit the Cox model of each process on `dt`, restricted to the rows at risk
# when the process has an at-risk condition. Returns a list named by process
# of list(fit, hazard), where `hazard` holds the positive baseline hazard
# increments `dhazard.<process>` (and the cumulative `hazard.<process>`) at
# each time.
.fit_cox_models <- function(
  dt,
  fit_types,
  at.risks,
  depend_time,
  recurrent,
  max_count,
  verbose
) {
  fits <- lapply(names(fit_types), function(name) {
    formula <- .cox_formula(
      fit_types[[name]][["model"]],
      depend_time,
      recurrent,
      max_count
    )
    data <- if (name %in% names(at.risks)) dt[at.risks[[name]](dt)] else dt
    fit <- coxph(
      as.formula(formula),
      data = data,
      control = coxph.control(timefix = FALSE)
    )
    if (verbose) {
      message("------------------------------")
      message(paste0("fit.type = ", name))
      .message_print(fit)
    }
    dhazard <- paste0("dhazard.", name)
    hazard <- suppressWarnings(setDT(basehaz(fit, centered = TRUE)))
    hazard[, (dhazard) := c(hazard[1], diff(hazard))]
    setnames(hazard, "hazard", paste0("hazard.", name))
    setcolorder(hazard, c("time", dhazard))
    list(fit = fit, hazard = hazard[hazard[[dhazard]] > 0])
  })
  stats::setNames(fits, names(fit_types))
}

# The model string of a process with the jump time variables `T.<name>` as
# factors and each recurrent process `z` replaced by its indicators
# `z.1 + ... + z.<max_count>`.
.cox_formula <- function(model, depend_time, recurrent, max_count) {
  for (T.var in sprintf("T.%s", names(depend_time))) {
    model <- gsub(T.var, paste0("as.factor(", T.var, ")"), model, fixed = TRUE)
  }
  for (process in recurrent) {
    model <- gsub(
      paste0("+", process),
      paste0("+", process, ".", seq_len(max_count), collapse = ""),
      model,
      fixed = TRUE
    )
  }
  model
}

# The hazard increments of the Cox fit `fit` on `newdata`: the baseline
# increments `dhazard` times the relative risk.
.cox_increment <- function(fit, newdata, dhazard) {
  as.numeric(dhazard * exp(predict(fit, newdata = newdata, type = "lp")))
}

# A Cox fit without the environment its formula was built in, which holds
# all of the fitting data and would be saved along with the fit. It keeps its
# model frame instead (only the rows and variables used), which is what
# survfit() and basehaz() need, e.g. in simevent::sim_model_from_fits().
.detach_cox_env <- function(fit) {
  fit$model <- stats::model.frame(fit)
  attr(attr(fit$model, "terms"), ".Environment") <- topenv()
  environment(fit$formula) <- topenv()
  attr(fit$terms, ".Environment") <- topenv()
  fit
}
