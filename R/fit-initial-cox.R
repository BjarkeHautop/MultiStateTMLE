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
  verbose
) {
  fits <- lapply(names(fit_types), function(name) {
    formula <- .cox_formula(
      fit_types[[name]][["model"]],
      depend_time
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

# The model string of a process with each term that is a jump time variable
# `T.<name>` made a factor, as each value of the variable is a state.
.cox_formula <- function(model, depend_time) {
  replacements <- list()
  for (T.var in sprintf("T.%s", names(depend_time))) {
    replacements[[T.var]] <- call("as.factor", as.name(T.var))
  }
  formula <- str2lang(model)
  formula[[3]] <- .replace_terms(formula[[3]], replacements)
  .deparse_one_line(formula)
}

# `expr`, a formula's right-hand side, with each term that is a name in
# `replacements` replaced by its value, looking into the formula operators
# that combine terms but not into other calls.
.replace_terms <- function(expr, replacements) {
  if (is.name(expr)) {
    replacement <- replacements[[as.character(expr)]]
    return(if (is.null(replacement)) expr else replacement)
  }
  operators <- c("+", "-", ":", "*", "/", "^", "%in%", "(")
  if (is.call(expr) && as.character(expr[[1]])[1] %in% operators) {
    for (i in seq_along(expr)[-1]) {
      expr[[i]] <- .replace_terms(expr[[i]], replacements)
    }
  }
  expr
}

# The hazard increments of the Cox fit `fit` on the data.table `newdata`:
# the baseline increments `dhazard` times the relative risk, in the rows
# `at.risk` (recycled; NULL for all), and 0 in the others. The rows not at
# risk are not predicted, as they can have covariate values the fit never
# saw, such as a factor level that only occurs where the process is not at
# risk.
.cox_increment <- function(fit, newdata, dhazard, at.risk = NULL) {
  if (is.null(at.risk)) {
    at.risk <- TRUE
  }
  at.risk <- rep_len(as.logical(at.risk), length(dhazard))
  if (anyNA(at.risk)) {
    stop("Internal error: missing at-risk indicators.", call. = FALSE)
  }
  out <- numeric(length(dhazard))
  rows <- which(at.risk)
  if (length(rows) > 0) {
    lp <- predict(fit, newdata = newdata[rows, ], type = "lp")
    out[rows] <- dhazard[rows] * exp(lp)
  }
  out
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
