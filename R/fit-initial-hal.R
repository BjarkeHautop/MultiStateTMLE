# Experimental HAL fits for fit_initial(). This code is untested: the HAL
# fits need a fit.hal() function on the search path, which MultiStateTMLE does
# not provide. It is kept apart from the Cox code path in fit-initial.R.

# Fit the HAL models of the processes with fit = "hal" on the expanded data
# `tmp.long`, returning the updated `tmp.long`, the fits and the state
# variables with non-zero HAL coefficients.
.hal_fit <- function(
  tmp.long,
  fit_types,
  which.hal,
  at.risks,
  hal_control,
  state.names,
  varnames,
  unique.times,
  A0.name,
  n,
  tau,
  max_count,
  cores,
  verbose
) {
  process.names <- names(fit_types)
  # Unpack the HAL settings into the names used by the HAL code below.
  cut.time <- hal_control$cut_time
  cut.one.way <- hal_control$cut_one_way
  cut.Tk <- hal_control$cut_tk
  cut.Tk.values <- hal_control$cut_tk_values
  two.way <- hal_control$two_way
  reduce.NK <- hal_control$reduce_nk
  hal.sl <- hal_control$sl
  if (is.null(hal.sl)) {
    hal.sl <- list(c(
      cut.one.way = cut.one.way,
      cut.Tk = cut.Tk,
      cut.Tk.values = cut.Tk.values,
      two.way = two.way
    ))
  }
  screen.two.way <- hal_control$screen_two_way
  lambda.cvs <- hal_control$lambda_cvs
  event.dependent.cv <- hal_control$event_dependent_cv
  npenalize.vars <- hal_control$npenalize_vars
  V <- hal_control$v_folds
  seed.hal <- hal_control$seed
  reduce.seed.dependence <- hal_control$reduce_seed_dependence
  penalize.time <- hal_control$penalize_time
  verbose.hal <- hal_control$verbose
  cv.glmnet <- FALSE
  max.count <- max_count
  use.cores <- cores

  deltas.hal <- vapply(
    fit_types[which.hal],
    function(fit.type) fit.type[["delta"]],
    numeric(1)
  )

  tmp.long[, idN := 1:.N, by = "id"]

  N.vars <- state.names
  T.vars <- grep("T.", varnames, value = TRUE)

  for (N.var in N.vars[N.vars %in% process.names]) {
    N.var.count <- tmp.long[, (get(N.var))[.N], by = "id"]
    NK <- ceiling(quantile(N.var.count[[2]], probs = 1 - reduce.NK))
    for (jnk in 1:NK) {
      T.vars <- unique(c(T.vars, paste0("T.", N.var, ".", jnk)))
    }
  }

  N.vars <- unique(c(N.vars, T.vars))

  baseline.vars <- setdiff(varnames, N.vars)

  ##T.vars2 <- unique(c(T.vars, unlist(lapply(hal.sl, function(hsl) hsl[grep("T.vars", names(hsl))]))))
  ##T.vars2 <- T.vars2[T.vars2 != ""]

  for (T.var in T.vars) {
    time.value <- as.numeric(strsplit(T.var, "\\.")[[1]][3])
    if (is.na(time.value)) {
      time.value <- 1
    }
    time.var <- strsplit(T.var, "\\.")[[1]][2]
    suppressWarnings(
      tmp.long[,
        (T.var) := min(time[
          c(get(time.var)[-1], get(time.var)[.N]) == time.value
        ]),
        by = "id"
      ]
    )
    tmp.long[tstart < get(T.var), (T.var) := 0]
  }

  tmp.long[,
    (varnames) := lapply(.SD, function(x) {
      if (is.character(x)) {
        return(as.numeric(as.factor(x)))
      } else if (!is.numeric(x)) {
        return(as.numeric(x))
      } else {
        return(x)
      }
    }),
    .SDcols = varnames
  ]

  grid.times <- unique.times[floor(seq(
    1,
    length(unique.times),
    length = cut.time + 2
  ))[-c(1, cut.time + 2)]]
  #grid.times <- grid.times[grid.times != 0.00001]

  tmp.long[,
    grid.period := as.numeric(cut(
      tstart,
      unique(c(0, grid.times, Inf)),
      include.lowest = TRUE,
      right = FALSE
    ))
  ]
  tmp.long[, grid.time := c(0, grid.times, Inf)[grid.period]]

  by.vars <- c("id", "grid.period")
  if (length(N.vars) > 0) {
    by.vars <- c(by.vars, N.vars)
  }

  for (delta.value in deltas.hal) {
    name.hal <- names(which.hal)[which.hal][deltas.hal == delta.value]
    if (name.hal %in% names(at.risks)) {
      tmp.long[,
        (paste0("at.risk.", delta.value)) := at.risks[[name.hal]](tmp.long)
      ]
    } else {
      tmp.long[, (paste0("at.risk.", delta.value)) := (time <= final.time)]
    }

    if (delta.value == 0) {
      tmp.long[,
        (paste0("ddd.", delta.value)) := sum(
          get(paste0("at.risk.", delta.value)) &
            time == final.time &
            delta == delta.value
        ),
        by = by.vars
      ]
    } else {
      tmp.long[,
        (paste0("ddd.", delta.value)) := sum(
          get(paste0("at.risk.", delta.value)) &
            (time <= final.time) * (delta == delta.value)
        ),
        by = by.vars
      ]
    }

    #tmp.long[, (paste0("risk.time.", delta.value)) := 0]
    tmp.long[
      time <= final.time & get(paste0("at.risk.", delta.value)),
      (paste0("risk.time.", delta.value)) := sum(time - tstart),
      by = by.vars
    ]
  }

  tmp.long.reduced <- unique(tmp.long[time <= final.time], by = by.vars)

  make.hal.formula <- function(
    cut.one.way = cut.one.way,
    cut.Tk = cut.Tk,
    cut.Tk.values = cut.Tk.values,
    two.way = two.way,
    N.vars = N.vars,
    T.vars = T.vars,
    max.count = max.count
  ) {
    grid.baseline.vars <- lapply(baseline.vars, function(var) {
      sort(unique(tmp.long[idN == 1, var, with = FALSE][floor(seq(
        1,
        n,
        length = cut.one.way
      ))])[[1]])
    })
    names(grid.baseline.vars) <- baseline.vars

    if (!all(T.vars == "")) {
      if (length(cut.Tk.values) > 0 && length(cut.Tk.values[[1]]) > 0) {
        grid.T.vars <- lapply(T.vars, function(T.var) {
          cut.Tk.values
        })
      } else {
        grid.T.vars <- lapply(T.vars, function(T.var) {
          tmp.var <- strsplit(gsub("T.", "", T.var), "\\.")[[1]]
          N.var <- tmp.var[1]
          N.count <- tmp.var[2]
          unique.T.var <- sort(unique(
            tmp.long[time <= tau & get(N.var) >= N.count][[T.var]][
              c(
                0,
                diff(tmp.long[time <= tau & get(N.var) >= N.count][[T.var]])
              ) >
                0
            ]
          ))
          return(unique.T.var[floor(seq(
            1,
            length(unique.T.var),
            length = min(cut.Tk, length(unique.T.var) / 3)
          ))])
        })
      }
      names(grid.T.vars) <- T.vars
    } else {
      T.vars <- NULL
    }

    if (all(N.vars == "")) {
      N.vars <- NULL
    }

    N.vars <- unique(c(
      setdiff(N.vars, grep("T.", N.vars, value = TRUE)),
      T.vars
    ))

    hal.formula.main <- paste0(
      "delta ~ -1 + ",
      ifelse(length(A0.name) > 0, paste0(A0.name, " + "), ""),
      paste0("(grid.time >=", grid.times, ")", collapse = "+"),
      " + ",
      paste0(
        sapply(setdiff(baseline.vars, A0.name), function(var) {
          paste0(
            "(",
            var,
            ">=",
            grid.baseline.vars[[var]],
            ")",
            collapse = "+"
          )
        }),
        collapse = "+"
      ),
      ifelse(
        length(setdiff(N.vars, T.vars)) > 0,
        paste0(
          " + ",
          paste0(
            sapply(setdiff(N.vars, T.vars), function(var) {
              paste0("(", var, ">=", 1:max.count, ")", collapse = "+")
            }),
            collapse = "+"
          )
        ),
        ""
      ),
      ifelse(
        length(T.vars) > 0,
        paste0(
          " + ",
          paste0(
            sapply(T.vars, function(var) {
              paste0("(", var, ">=", grid.T.vars[[var]], ")", collapse = "+")
            }),
            collapse = "+"
          )
        ),
        ""
      )
    )

    if (length(two.way) > 0 && length(two.way[[1]]) > 0) {
      two.way.matrix <- do.call(
        "rbind",
        lapply(two.way, function(two.way.jj) {
          two.way.inner <- unique(
            do.call(
              "rbind",
              combn(
                unlist(lapply(seq_along(two.way.jj[1:2]), function(ii) {
                  two.way.var <- two.way.jj[ii]
                  k1 <- as.numeric(two.way.jj[3])
                  k2 <- ifelse(
                    length(two.way.jj) >= 4,
                    as.numeric(two.way.jj[4]),
                    k1
                  )
                  kk <- if (ii == 1) k1 else k2
                  if (two.way.var == "time") {
                    return(paste0(
                      "(grid.time >=",
                      grid.times[seq(1, length(grid.times), length = kk)],
                      ")"
                    ))
                  } else if (two.way.var %in% baseline.vars) {
                    return(paste0(
                      "(",
                      two.way.var,
                      ">=",
                      grid.baseline.vars[[two.way.var]][seq(
                        1,
                        length(grid.baseline.vars[[two.way.var]]),
                        length = min(
                          length(grid.baseline.vars[[two.way.var]]),
                          kk
                        )
                      )],
                      ")"
                    ))
                  } else if (two.way.var %in% setdiff(N.vars, T.vars)) {
                    return(paste0(
                      "(",
                      two.way.var,
                      ">=",
                      (1:max.count)[seq(1, max.count, length = kk)],
                      ")"
                    ))
                  } else if (two.way.var %in% T.vars) {
                    return(paste0(
                      "(",
                      two.way.var,
                      ">=",
                      grid.T.vars[[two.way.var]][seq(
                        1,
                        length(grid.T.vars[[two.way.var]]),
                        length = min(length(grid.T.vars[[two.way.var]]), kk)
                      )],
                      ")"
                    ))
                  } else if (two.way.var == A0.name) {
                    return(A0.name)
                  }
                })),
                2,
                simplify = FALSE
              )
            )
          )
          two.way.inner <-
            two.way.inner[
              sapply(
                strsplit(
                  as.character(two.way.inner[, 1]),
                  ">=",
                  fixed = TRUE
                ),
                function(twi) {
                  gsub("\\(", "", twi[1])
                }
              ) !=
                sapply(
                  strsplit(
                    as.character(two.way.inner[, 2]),
                    ">=",
                    fixed = TRUE
                  ),
                  function(twi) {
                    gsub("\\(", "", twi[1])
                  }
                ),
            ]
          if (
            startsWith(two.way.jj[1], "T.") &&
              startsWith(two.way.jj[2], "T.")
          ) {
            two.way.inner <-
              two.way.inner[
                sapply(
                  strsplit(
                    as.character(two.way.inner[, 1]),
                    ">=",
                    fixed = TRUE
                  ),
                  function(twi) {
                    as.numeric(gsub(")", "", twi[2]))
                  }
                ) <
                  sapply(
                    strsplit(
                      as.character(two.way.inner[, 2]),
                      ">=",
                      fixed = TRUE
                    ),
                    function(twi) {
                      as.numeric(gsub(")", "", twi[2]))
                    }
                  ),
              ]
          }
          if (
            two.way.jj[1] == "time" &&
              startsWith(two.way.jj[2], "T.")
          ) {
            two.way.inner <-
              two.way.inner[
                sapply(
                  strsplit(
                    as.character(two.way.inner[, 1]),
                    ">=",
                    fixed = TRUE
                  ),
                  function(twi) {
                    as.numeric(gsub(")", "", twi[2]))
                  }
                ) >
                  sapply(
                    strsplit(
                      as.character(two.way.inner[, 2]),
                      ">=",
                      fixed = TRUE
                    ),
                    function(twi) {
                      as.numeric(gsub(")", "", twi[2]))
                    }
                  ),
              ]
          }
          if (
            two.way.jj[2] == "time" &&
              startsWith(two.way.jj[1], "T.")
          ) {
            two.way.inner <-
              two.way.inner[
                sapply(
                  strsplit(
                    as.character(two.way.inner[, 1]),
                    ">=",
                    fixed = TRUE
                  ),
                  function(twi) {
                    as.numeric(gsub(")", "", twi[2]))
                  }
                ) <
                  sapply(
                    strsplit(
                      as.character(two.way.inner[, 2]),
                      ">=",
                      fixed = TRUE
                    ),
                    function(twi) {
                      as.numeric(gsub(")", "", twi[2]))
                    }
                  ),
              ]
          }
          return(two.way.inner[
            apply(two.way.inner, 1, function(two.way.inner.row) {
              !((length(grep(two.way.jj[1], two.way.inner.row[1])) > 0 &
                length(grep(two.way.jj[1], two.way.inner.row[2])) > 0) |
                (length(grep(two.way.jj[2], two.way.inner.row[1])) > 0 &
                  length(grep(two.way.jj[2], two.way.inner.row[2])) > 0))
            }),
          ])
        })
      )
      hal.formula.main <- paste0(
        hal.formula.main,
        " + ",
        paste0(two.way.matrix[, 1], ":", two.way.matrix[, 2], collapse = "+")
      )
    }
    return(hal.formula.main)
  }

  make.hal <- function(hal.jj) {
    if (
      sum(
        names(hal.jj) %in%
          c(
            "cut.one.way",
            "max.count",
            "T.vars",
            "N.vars",
            "cut.Tk",
            "cut.Tk.values",
            "two.way"
          )
      ) >
        0
    ) {
      if (length(grep("max.count", hal.jj)) > 0) {
        hal.jj["max.count"] <- min(hal.jj["max.count"], NK)
      }
      return(hal.jj)
    }

    if (sum(names(hal.jj) %in% c("components", "interactions")) == 0) {
      split.hal.jj <- strsplit(hal.jj, "_")
      if (length(split.hal.jj[[1]]) > 1) {
        use.times <- gsub("time", "", split.hal.jj[[1]][2])
      } else {
        use.times <- NULL
      }

      hal.jj <- split.hal.jj[[1]][1]

      if (hal.jj == "default") {
        hal.jj <- list(
          components = list(
            baseline = list(complexity = "medium"),
            history = list(complexity = "medium", time.mode = "auto")
          ),
          interactions = list(
            time.with = c("history", A0.name),
            history.with = c("history"),
            baseline.with = c("baseline", A0.name)
          )
        )
      } else if (hal.jj == "history.interactions") {
        hal.jj <- list(
          components = list(
            baseline = list(complexity = "medium"),
            history = list(complexity = "medium", time.mode = "auto")
          ),
          interactions = list(
            time.with = c("history"),
            history.with = c("history")
          )
        )
      } else if (hal.jj == "dense") {
        hal.jj <- list(
          components = list(
            baseline = list(complexity = "high"),
            history = list(complexity = "high", time.mode = "auto")
          ),
          interactions = list(
            time.with = c("history", A0.name),
            baseline.with = c("baseline", A0.name)
          )
        )
      } else if (hal.jj == "full2") {
        hal.jj <- list(
          components = list(
            baseline = list(complexity = "high"),
            history = list(complexity = "high", time.mode = "auto")
          ),
          interactions = list(
            time.with = c("history", A0.name, "baseline"),
            history.with = c("history", "baseline", A0.name),
            baseline.with = c("baseline", A0.name)
          )
        )
      } else if (hal.jj == "full") {
        hal.jj <- list(
          components = list(
            baseline = list(complexity = "high"),
            history = list(complexity = "high", time.mode = "auto")
          ),
          interactions = list(
            time.with = c("history", A0.name),
            history.with = c("history"),
            baseline.with = c("baseline", A0.name)
          )
        )
      } else if (hal.jj == "no.interactions") {
        hal.jj <- list(
          components = list(
            baseline = list(complexity = "medium"),
            history = list(complexity = "medium", time.mode = "auto")
          )
        )
      } else if (hal.jj == "baseline") {
        hal.jj <- list(
          components = list(
            baseline = list(complexity = "medium"),
            history = list(complexity = "none", time.mode = "none")
          )
        )
      } else if (hal.jj == "sparse") {
        hal.jj <- list(
          components = list(
            baseline = list(complexity = "medium"),
            history = list(complexity = "medium", time.mode = "none")
          )
        )
      } else {
        message("NB: the choice of hal not supported - sparse will be fitted")
        hal.jj <- list(
          components = list(
            baseline = list(complexity = "medium"),
            history = list(complexity = "medium", time.mode = "none")
          )
        )
      }

      if (length(use.times) > 0) {
        hal.jj[["components"]][length(hal.jj[["components"]]) + 1] <-
          list(list(use.times = use.times))
        names(hal.jj[["components"]])[length(hal.jj[["components"]])] <-
          "time"
      }
    }

    if (sum(names(hal.jj) %in% c("components", "interactions")) > 0) {
      hal.out <- c()
      baseline <- hal.jj[["components"]][["baseline"]]
      if (length(baseline) > 0) {
        if (baseline$complexity == "medium") {
          hal.out <- c(hal.out, cut.one.way = 20)
        } else if (baseline$complexity == "high") {
          hal.out <- c(hal.out, cut.one.way = 35)
        } else {
          hal.out <- c(hal.out, cut.one.way = 10)
        }
      }
      history <- hal.jj[["components"]][["history"]]
      if (length(history) > 0) {
        if (history$complexity == "medium") {
          hal.out <- c(hal.out, max.count = min(3, NK), cut.Tk = 5)
        } else if (history$complexity == "high") {
          hal.out <- c(hal.out, max.count = min(10, NK), cut.Tk = 10)
        } else if (history$complexity == "low") {
          hal.out <- c(hal.out, max.count = min(2, NK), cut.Tk = 2)
        } else {
          hal.out <- c(hal.out, max.count = 0, cut.Tk = 0)
        }
        if (!("time.mode" %in% names(history))) {
          history$time.mode <- "auto"
        }
        if (history$time.mode == "auto") {
          hal.out <- c(
            hal.out,
            T.vars = grep(
              paste0(".", 1:hal.out["max.count"], collapse = "|"),
              T.vars,
              value = TRUE
            )
          )
        } else if (history$time.mode == "first" && hal.out["max.count"] >= 1) {
          hal.out <- c(
            hal.out,
            T.vars = grep(
              grep("\\.1$", T.vars, value = TRUE),
              T.vars,
              value = TRUE
            )
          )
        } else {
          hal.out <- c(hal.out, T.vars = "")
        }
      }
      ## this is test, not recommended
      time <- hal.jj[["components"]][["time"]]
      if (length(time) > 0) {
        if ("complexity" %in% names(time)) {
          if (time$complexity == "medium") {
            hal.out <- c(hal.out, penalize.time = 0.4)
          } else if (time$complexity == "high") {
            hal.out <- c(hal.out, penalize.time = 0)
          } else if (time$complexity == "none") {
            hal.out <- c(hal.out, penalize.time = 1)
          } else {
            hal.out <- c(hal.out, penalize.time = 0.7)
          }
        }
        if ("use.times" %in% names(time)) {
          hal.out <- c(hal.out, use.times = time$use.times)
        }
      }

      if (sum(names(hal.jj) %in% c("interactions")) > 0) {
        two.way.out <- list()

        time.with <- hal.jj[["interactions"]][["time.with"]]
        if (length(time.with) > 0) {
          if ("history" %in% time.with) {
            for (T.var in hal.out[grep("T.vars", names(hal.out))]) {
              two.way.out[[length(two.way.out) + 1]] <-
                c("time", T.var, 30, 10)
            }
            for (N.var in setdiff(N.vars, grep("T.", N.vars, value = TRUE))) {
              two.way.out[[length(two.way.out) + 1]] <-
                c("time", N.var, 10, 10)
            }
          }
          if (length(A0.name) > 0) {
            if (A0.name %in% time.with) {
              two.way.out[[length(two.way.out) + 1]] <-
                c("time", A0.name, 10)
            }
          }
          if ("baseline" %in% time.with) {
            for (baseline.var in setdiff(baseline.vars, A0.name)) {
              two.way.out[[length(two.way.out) + 1]] <-
                c("time", baseline.var, 10)
            }
          }
        }

        ### fix: need to make this:
        history.with <- hal.jj[["interactions"]][["history.with"]]
        if (length(history.with) > 0) {
          if ("history" %in% history.with) {
            T.vars <- hal.out[grep("T.vars", names(hal.out))]
            length.T.vars <- length(T.vars)
            if (length.T.vars > 1) {
              for (tt1 in 1:(length.T.vars - 1)) {
                for (tt2 in (tt1 + 1):length(T.vars)) {
                  two.way.out[[length(two.way.out) + 1]] <-
                    c(T.vars[tt1], T.vars[tt2], min(10, hal.out["cut.Tk"]))
                }
              }
            }
          }
          if (length(A0.name) > 0) {
            if (A0.name %in% history.with) {
              for (T.var in hal.out[grep("T.vars", names(hal.out))]) {
                two.way.out[[length(two.way.out) + 1]] <-
                  c(A0.name, T.var, 10)
              }
              for (N.var in setdiff(
                N.vars,
                grep("T.", N.vars, value = TRUE)
              )) {
                two.way.out[[length(two.way.out) + 1]] <-
                  c(A0.name, N.var, 10)
              }
            }
          }
          if ("baseline" %in% history.with) {
            for (baseline.var in setdiff(baseline.vars, A0.name)) {
              for (T.var in hal.out[grep("T.vars", names(hal.out))]) {
                two.way.out[[length(two.way.out) + 1]] <-
                  c(baseline.var, T.var, 10)
              }
              for (N.var in setdiff(
                N.vars,
                grep("T.", N.vars, value = TRUE)
              )) {
                two.way.out[[length(two.way.out) + 1]] <-
                  c(baseline.var, N.var, 10)
              }
            }
          }
        }

        baseline.with <- hal.jj[["interactions"]][["baseline.with"]]
        if (length(baseline.with) > 0) {
          if ("baseline" %in% baseline.with) {
            baseline.covars <- setdiff(baseline.vars, A0.name)
            length.baseline <- length(baseline.covars)
            if (length.baseline > 1) {
              for (baseline.j1 in 1:(length.baseline - 1)) {
                for (baseline.j2 in (baseline.j1 + 1):length(
                  baseline.covars
                )) {
                  two.way.out[[length(two.way.out) + 1]] <-
                    c(
                      baseline.covars[baseline.j1],
                      baseline.covars[baseline.j2],
                      10
                    )
                }
              }
            }
          }
          if (length(A0.name) > 0) {
            if (A0.name %in% baseline.with) {
              baseline.covars <- setdiff(baseline.vars, A0.name)
              for (baseline.j1 in seq_along(baseline.covars)) {
                two.way.out[[length(two.way.out) + 1]] <-
                  c(baseline.covars[baseline.j1], A0.name, 10)
              }
            }
          }
        }
      } else {
        two.way.out <- list()
      }

      hal.out <- c(hal.out, two.way = unique(two.way.out))
    }

    return(hal.out)
  }

  ## FIX! need to make this work
  fit.hals.list <- lapply(seq_along(hal.sl), function(hal.jj) {
    hal.parameters <- make.hal(hal.jj = hal.sl[[hal.jj]])

    if (length(grep("penalize.time", names(hal.parameters))) > 0) {
      penalize.time <- as.numeric(hal.parameters[grep(
        "penalize.time",
        names(hal.parameters)
      )])
    }

    if (length(grep("use.times", names(hal.parameters))) > 0) {
      use.times <- as.numeric(hal.parameters[grep(
        "use.times",
        names(hal.parameters)
      )])
    } else {
      use.times <- 1
    }

    if (verbose) {
      message(paste0(hal.jj, "/", length(hal.sl), " (hal-sl)"))
    }

    hal.formula.main <- make.hal.formula(
      cut.one.way = ifelse(
        "cut.one.way" %in% names(hal.parameters),
        as.numeric(hal.parameters[["cut.one.way"]]),
        cut.one.way
      ),
      cut.Tk = if (length(grep("cut.Tk", names(hal.parameters))) > 0) {
        as.numeric(hal.parameters[grep("cut.Tk", names(hal.parameters))])
      } else {
        cut.Tk
      },
      cut.Tk.values = if (
        length(grep("cut.Tk.values", names(hal.parameters))) > 0
      ) {
        hal.parameters[grep("cut.Tk.values", names(hal.parameters))]
      } else {
        cut.Tk.values
      },
      two.way = if (any(grepl("two.way", names(hal.parameters)))) {
        which.two.way <- grep("two.way", names(hal.parameters))
        lapply(which.two.way, function(tw) hal.parameters[[tw]])
      } else {
        two.way
      },
      N.vars = if (length(grep("N.vars", names(hal.parameters))) > 0) {
        hal.parameters[grep("N.vars", names(hal.parameters))]
      } else {
        N.vars
      },
      T.vars = if (length(grep("T.vars", names(hal.parameters))) > 0) {
        hal.parameters[grep("T.vars", names(hal.parameters))]
      } else {
        T.vars
      },
      max.count = ifelse(
        "max.count" %in% names(hal.parameters),
        as.numeric(hal.parameters[["max.count"]]),
        max.count
      )
    )

    X <- Matrix(
      model.matrix(formula(hal.formula.main), data = tmp.long.reduced),
      sparse = TRUE
    )

    X <- X[, !(colnames(X) %in% grep("FALSE", colnames(X), value = TRUE))]

    ##col_ones <- Matrix::colMeans(X != 0)
    ##X <- X[, col_ones >= min.no.of.ones]

    x.vector <- hash_sparse_rows_dgC(X)
    tmp.long.reduced[, x := x.vector]

    if (length(seed.hal) == 0) {
      seed.hal <- sample(1e8, 1)
    }

    fit.hals <- lapply(deltas.hal, function(delta.value) {
      ##message("-------")
      if (verbose) {
        message(paste0("delta = ", delta.value))
      }
      ##message("-------")

      tmp.long.reduced[, risk.time := get(paste0("risk.time.", delta.value))]

      first.hal <- fit.hal(
        hal.pseudo.dt = tmp.long.reduced[tmp.long.reduced[[paste0(
          "at.risk.",
          delta.value
        )]]],
        X.hal = X[tmp.long.reduced[[paste0("at.risk.", delta.value)]], ],
        delta.var = "delta",
        delta.value = delta.value,
        time.var = "time",
        treatment = A0.name,
        lambda.cvs = lambda.cvs,
        penalize.time = if (penalize.time < 1) {
          penalize.time
        } else {
          screen.two.way | penalize.time
        }, #FALSE,
        penalize.treatment = screen.two.way, #FALSE,
        npenalize.vars = npenalize.vars,
        event.dependent.cv = event.dependent.cv,
        reduce.seed.dependence = reduce.seed.dependence,
        return.train.fit = FALSE,
        use.times = use.times,
        check.support.for.basis.functions = 0.01,
        V = V,
        parallelize.cve = min(detectCores() - 1, use.cores),
        cv.glmnet = cv.glmnet,
        verbose = verbose.hal,
        return.cve = TRUE,
        seed = seed.hal
      )

      if (screen.two.way) {
        get.coef <- coef(first.hal$hal.fit, s = first.hal$lambda)
        coef.dt <- data.table(
          non.zero = get.coef@Dimnames[[1]][get.coef@i + 1],
          coef = get.coef@x
        ) #[non.zero %in% c("Aobs", "Y.dummy >= 1TRUE")]
        first.hal.coefs <- coef.dt[["coef"]]
        which.two.way <- names(first.hal.coefs) <- gsub(
          "TRUE|FALSE",
          "",
          paste0(coef.dt[["non.zero"]])
        )
        two.way <- do.call(
          "rbind",
          combn(
            which.two.way[which.two.way != "(Intercept)"],
            2,
            simplify = FALSE
          )
        )
        two.way <- two.way[
          apply(two.way, 1, function(two.way.row) {
            (strsplit(two.way.row[1], ">=")[[1]][1] !=
              strsplit(two.way.row[2], ">=")[[1]][1])
          }),
        ]
        hal.formula.main.2 <-
          paste0(
            hal.formula.main,
            " + ",
            paste0(
              "(",
              two.way[, 1],
              "):(",
              two.way[, 2],
              ")",
              collapse = "+"
            )
          )
        X.2 <- Matrix(
          model.matrix(formula(hal.formula.main.2), data = tmp.long.reduced),
          sparse = TRUE
        )
        return(
          fit.hal(
            hal.pseudo.dt = tmp.long.reduced[tmp.long.reduced[[paste0(
              "at.risk.",
              delta.value
            )]]],
            X.hal = X.2[tmp.long.reduced[[paste0("at.risk.", delta.value)]]],
            delta.var = "delta",
            delta.value = delta.value,
            time.var = "time",
            treatment = A0.name,
            lambda.cvs = lambda.cvs,
            penalize.time = FALSE | penalize.time,
            penalize.treatment = FALSE,
            event.dependent.cv = event.dependent.cv,
            reduce.seed.dependence = reduce.seed.dependence,
            return.train.fit = FALSE,
            use.times = use.times,
            check.support.for.basis.functions = TRUE,
            V = V,
            parallelize.cve = min(detectCores() - 1, use.cores),
            cv.glmnet = cv.glmnet,
            verbose = verbose.hal,
            return.cve = TRUE,
            seed = seed.hal
          )
        )
      } else {
        return(first.hal)
      }
    })
  })

  if (length(fit.hals.list) == 1) {
    fit.hals <- fit.hals.list[[1]]
    if (!verbose.hal && verbose) {
      lapply(fit.hals, function(fh) {
        message("--------------------------------------------")
        message(paste0("delta = ", fh$delta.value))
        message("--------------------------------------------")
        .message_print(coef(fh$hal.fit, s = fh$lambda.cv))
      })
    }
  } else if (length(fit.hals.list) > 1) {
    fit.hals <- lapply(
      seq_along(fit.hals.list[[1]]),
      function(fh.process.index) {
        if (verbose) {
          message("--------------------------------------------")
          message(paste0(
            "delta = ",
            fit.hals.list[[1]][[fh.process.index]]$delta.value
          ))
          message("--------------------------------------------")
        }
        # here:
        inner.cve <- sapply(fit.hals.list, function(fh) {
          fh[[fh.process.index]]$cve
        })
        picked <- (seq_along(inner.cve))[inner.cve == min(inner.cve)]
        if (verbose) {
          .message_print(inner.cve)
        }
        if (verbose) {
          message(paste0("picked hal: ", picked))
        }
        fh <- fit.hals.list[[picked]][[fh.process.index]]
        if (verbose) {
          .message_print(coef(fh$hal.fit, s = fh$lambda.cv))
        }
        return(fh)
      }
    )
  }

  tmp.long[, risk.time := time - tstart]

  hal.coefs <- lapply(fit.hals, function(fh) {
    delta <- fh$delta.value
    get.coef <- coef(fh$hal.fit, s = fh$lambda)
    coef.dt <- data.table(
      non.zero = get.coef@Dimnames[[1]][get.coef@i + 1],
      coef = get.coef@x
    ) #[non.zero %in% c("Aobs", "Y.dummy >= 1TRUE")]
    coefs <- coef.dt[["coef"]]
    names(coefs) <- paste0(coef.dt[["non.zero"]])
    return(coefs)
  })

  non.zero.vars <- unique(unlist(sapply(
    hal.coefs[deltas.hal != 0],
    function(hal.coef) names(hal.coef)
  )))

  non.zero.state.vars <- lapply(N.vars, function(N.var) {
    tmp.vars <- gsub("TRUE", "", grep(N.var, non.zero.vars, value = TRUE))
    tmp.interaction <- grep(":", tmp.vars, value = TRUE)
    if (length(tmp.interaction) > 0) {
      tmp.interaction.vars <- unlist(lapply(
        strsplit(tmp.interaction, ":"),
        function(interaction.split) {
          grep(N.var, interaction.split, value = TRUE)
        }
      ))
      tmp.vars <- unique(c(
        setdiff(tmp.vars, tmp.interaction),
        tmp.interaction.vars
      ))
    }
    return(sort(unique(tmp.vars[
      !tmp.vars %in% grep(paste0("T.", N.var), tmp.vars, value = TRUE)
    ])))
  })
  names(non.zero.state.vars) <- N.vars
  ## the following is to make sure that N.x >= y is always included when T.x.y is
  for (T.var in T.vars) {
    if (length(non.zero.state.vars[[T.var]]) > 0) {
      time.value <- as.numeric(strsplit(T.var, "\\.")[[1]][3])
      if (is.na(time.value)) {
        time.value <- 1
      }
      time.var <- strsplit(T.var, "\\.")[[1]][2]
      if (
        all(
          non.zero.state.vars[[time.var]] !=
            paste0(time.var, " >= ", time.value)
        )
      ) {
        non.zero.state.vars[[time.var]] <-
          sort(c(
            non.zero.state.vars[[time.var]],
            paste0(time.var, " >= ", time.value)
          ))
      }
    }
  }
  ##message("FIX: CHECK this below")
  non.zero.state.vars <- non.zero.state.vars[(seq_along(
    non.zero.state.vars
  ))[sapply(non.zero.state.vars, function(non.zero) {
    length(non.zero) > 0
  })]]

  list(
    tmp.long = tmp.long,
    fits = fit.hals,
    which.hal = which.hal,
    N.vars = N.vars,
    T.vars = T.vars,
    non.zero.state.vars = non.zero.state.vars
  )
}

# The HAL design matrices on `tmp.long`, split into the columns that do not
# depend on the state or treatment (static), those that depend on the state
# (dynamic), on the treatment (A0) or on both (A0.dynamic).
.hal_design <- function(hal, tmp.long, A0.name) {
  fit.hals <- hal$fits
  N.vars <- hal$N.vars

  hal.vars.list <- lapply(seq_along(fit.hals), function(kk) {
    test.get.beta <- try(
      out.beta <- fit.hals[[kk]][["hal.fit"]]$beta@Dimnames[[1]],
      silent = TRUE
    )
    if (any(class(test.get.beta) == "try-error")) {
      out.beta <- fit.hals[[kk]][["hal.fit"]]$glmnet$beta@Dimnames[[1]]
    }
    return(out.beta)
  })

  hal.vars <- unique(c(unlist(hal.vars.list)))

  if (length(A0.name) > 0 && any(grepl(":", hal.vars))) {
    hal.vars.A0.dynamic <- unlist(sapply(N.vars, function(N.var) {
      tmp.vars <- #gsub("TRUE", "", grep(N.var, unlist(hal.vars.list), value = TRUE))
        grep(N.var, unlist(hal.vars.list), value = TRUE)
      tmp.vars <- tmp.vars[
        tmp.vars %in% grep(A0.name, unlist(hal.vars.list), value = TRUE)
      ]
      return(unique(tmp.vars[
        !tmp.vars %in% grep(paste0("T.", N.var), tmp.vars, value = TRUE)
      ]))
    }))
  } else {
    hal.vars.A0.dynamic <- NULL
  }

  hal.vars.dynamic <- setdiff(
    unlist(sapply(N.vars, function(N.var) {
      tmp.vars <- #gsub("TRUE", "", grep(N.var, unlist(hal.vars.list), value = TRUE))
        grep(N.var, unlist(hal.vars.list), value = TRUE)
      return(unique(tmp.vars[
        !tmp.vars %in% grep(paste0("T.", N.var), tmp.vars, value = TRUE)
      ]))
    })),
    hal.vars.A0.dynamic
  )

  if (length(A0.name) > 0) {
    hal.vars.A0 <- setdiff(
      grep(
        A0.name,
        hal.vars[!(hal.vars %in% grep("FALSE", hal.vars, value = TRUE))],
        value = TRUE
      ),
      hal.vars.A0.dynamic
    )
  } else {
    hal.vars.A0 <- NULL
  }

  hal.vars.static <- setdiff(
    hal.vars,
    c(
      hal.vars.dynamic,
      hal.vars.A0,
      hal.vars.A0.dynamic,
      unique(unlist(sapply(N.vars, function(N.var) {
        grep(N.var, unlist(hal.vars.list), value = TRUE)
      })))
    )
  )

  hal$vars.list <- hal.vars.list
  hal$vars.static <- hal.vars.static
  hal$vars.dynamic <- hal.vars.dynamic
  hal$vars.A0 <- hal.vars.A0
  hal$vars.A0.dynamic <- hal.vars.A0.dynamic
  hal$X.static <- .hal_matrix(hal.vars.static, tmp.long)
  hal$X.dynamic <- .hal_matrix(hal.vars.dynamic, tmp.long)
  hal$X.A0 <- if (length(A0.name) > 0) {
    .hal_matrix(hal.vars.A0, tmp.long)
  } else {
    .hal_empty_matrix(tmp.long)
  }
  hal$X.A0.dynamic <- if (length(hal.vars.A0.dynamic) > 0) {
    .hal_matrix(hal.vars.A0.dynamic, tmp.long)
  } else {
    .hal_empty_matrix(tmp.long)
  }
  hal
}

# Rebuild the treatment parts of the HAL design after the treatment column of
# `data` has been set to an arm.
.hal_design_set_arm <- function(hal, data) {
  hal$X.A0 <- .hal_matrix(hal$vars.A0, data)
  if (length(hal$vars.A0.dynamic) > 0) {
    hal$X.A0.dynamic <- .hal_matrix(hal$vars.A0.dynamic, data)
  }
  hal
}

# Rebuild the state parts of the HAL design for the state columns of `data`.
.hal_design_set_state <- function(hal, data) {
  hal$X.dynamic <- .hal_matrix(hal$vars.dynamic, data)
  if (length(hal$vars.A0.dynamic) > 0) {
    hal$X.A0.dynamic <- .hal_matrix(hal$vars.A0.dynamic, data)
  }
  hal
}

# The HAL hazard increments of process `name` on the current design.
.hal_predict <- function(hal, name, risk.time) {
  k <- match(name, names(hal$which.hal)[hal$which.hal])
  fit <- hal$fits[[k]]
  X <- cbind(hal$X.static, hal$X.dynamic, hal$X.A0.dynamic, hal$X.A0)
  exp(predict(
    fit[["hal.fit"]],
    X[, hal$vars.list[[k]]],
    newoffset = 0,
    s = fit[["lambda.cv"]]
  )) *
    risk.time
}

# The sparse design matrix of the HAL basis functions `vars` on `data`.
.hal_matrix <- function(vars, data) {
  terms <- paste0("(", gsub("FALSE|TRUE", "", gsub(":", "):(", vars)), ")")
  X <- Matrix(
    model.matrix(
      formula(paste0("delta~-1+", paste(terms, collapse = "+"))),
      data = data
    ),
    sparse = TRUE
  )
  .flip_interactions(X, vars)
}

.hal_empty_matrix <- function(data) {
  Matrix(0, nrow = nrow(data), ncol = 0, sparse = TRUE)
}

# Rename the interaction columns of X.mat that model.matrix() named in the
# other order than var.vec.
.flip_interactions <- function(X.mat, var.vec) {
  cols.diff <- setdiff(var.vec, colnames(X.mat))
  if (length(cols.diff) > 0) {
    cols.diff.rev <- sapply(cols.diff, function(x) {
      if (!grepl(":", x)) {
        return(x)
      }
      parts <- strsplit(x, ":", fixed = TRUE)[[1]]
      paste(rev(parts), collapse = ":")
    })
    cols.diff.TF <- grep("TRUE|FALSE", cols.diff, value = TRUE)
    if (length(cols.diff.TF) > 0) {
      idx <- match(cols.diff.rev, colnames(X.mat))
      valid <- !is.na(idx)
      colnames(X.mat)[idx[valid]] <- cols.diff[valid]
    }
    if (length(cols.diff.TF) < length(cols.diff)) {
      idx <- match(cols.diff.rev, gsub("TRUE|FALSE", "", colnames(X.mat)))
      valid <- !is.na(idx)
      colnames(X.mat)[idx[valid]] <- cols.diff[valid]
    }
  }
  X.mat
}

# The state space from the state variables with non-zero HAL coefficients,
# or NULL if there are none. Also recodes the jump time columns of
# `tmp.long` to the HAL knots. Returns list(depend.matrix, tmp.long).
.hal_state_space <- function(hal, tmp.long) {
  non.zero.state.vars <- hal$non.zero.state.vars
  N.vars <- hal$N.vars
  T.vars <- hal$T.vars
  if (length(non.zero.state.vars) == 0) {
    return(list(depend.matrix = NULL, tmp.long = tmp.long))
  }
  grid.list <- lapply(
    seq_along(non.zero.state.vars),
    function(non.zero.jj) {
      non.zero.var <- names(non.zero.state.vars)[non.zero.jj]
      if (non.zero.var %in% setdiff(N.vars, T.vars)) {
        return(
          0:max(as.numeric(gsub(
            paste0(non.zero.var, ">="),
            "",
            gsub(" ", "", non.zero.state.vars[[non.zero.jj]])
          )))
        )
      } else {
        return(c(
          0,
          sapply(
            seq_along(non.zero.state.vars[[non.zero.var]]),
            function(sss) {
              return(as.numeric(gsub(
                paste0(paste0(non.zero.var, ">=")),
                "",
                gsub(" ", "", non.zero.state.vars[[non.zero.var]][sss])
              )))
            }
          )
        ))
      }
    }
  )
  depend.matrix <- as.data.table(do.call(expand.grid, grid.list))
  setnames(depend.matrix, names(non.zero.state.vars))
  #names(depend.matrix) <- names(non.zero.state.vars)
  T.vars <- grep("T.", names(non.zero.state.vars), value = TRUE)
  if (length(T.vars) > 0) {
    for (T.var in T.vars) {
      time.var <- strsplit(T.var, "\\.")[[1]][2]
      time.value <- as.numeric(strsplit(T.var, "\\.")[[1]][3])
      if (is.na(time.value)) {
        time.value <- 1
      }
      depend.matrix <- depend.matrix[
        !(depend.matrix[[time.var]] <= time.value - 1 &
          depend.matrix[[T.var]] > 0)
      ]
    }
  }

  for (non.zero.jj in (seq_along(non.zero.state.vars))[
    names(non.zero.state.vars) %in% T.vars
  ]) {
    non.zero.var <- names(non.zero.state.vars)[non.zero.jj]
    non.zero.var.states <- non.zero.state.vars[[non.zero.var]]
    tmp.long[
      eval(parse(text = paste0("!(", non.zero.var.states[1], ")"))),
      (paste0(non.zero.var, "_1")) := 0
    ]
    tmp.long[
      eval(parse(
        text = paste0(
          "!(",
          gsub(non.zero.var, "time", non.zero.var.states[1]),
          ")"
        )
      )),
      (gsub("T.", "time.", non.zero.var)) := 0
    ]
    for (each.jj in seq_along(non.zero.var.states)) {
      value.jj <- as.numeric(gsub(
        paste0(paste0(non.zero.var, ">=")),
        "",
        gsub(" ", "", non.zero.state.vars[[non.zero.var]][each.jj])
      ))
      tmp.long[
        eval(parse(text = non.zero.var.states[each.jj])),
        (paste0(non.zero.var, "_1")) := value.jj
      ]
      tmp.long[
        eval(parse(
          text = gsub(non.zero.var, "time", non.zero.var.states[each.jj])
        )),
        (gsub("T.", "time.", non.zero.var)) := value.jj
      ]
    }
    #tmp.long[[paste0(non.zero.var)]] <- tmp.long[[paste0(non.zero.var, "_1")]]
    #tmp.long[[paste0(non.zero.var, "_1")]] <- NULL
    tmp.long[, (non.zero.var) := get(paste0(non.zero.var, "_1"))]
    tmp.long[, (paste0(non.zero.var, "_1")) := NULL]
  }
  list(depend.matrix = depend.matrix, tmp.long = tmp.long)
}

# The jump time variables of the HAL state space.
.hal_time_vars <- function(hal) {
  grep("T.", names(hal$non.zero.state.vars), value = TRUE)
}

# Hash each row of a sparse matrix by its sparsity pattern and values, as a
# fast way to find unique rows.
hash_sparse_rows_dgC <- function(M) {
  stopifnot(inherits(M, "dgCMatrix"))

  p <- M@p
  i <- M@i
  x <- M@x
  nr <- M@Dim[1]

  # storage: list of integer/values per row
  idx_list <- vector("list", nr)
  val_list <- vector("list", nr)

  # loop over columns, add entries to each row
  for (col in seq_len(M@Dim[2])) {
    start <- p[col] + 1
    end <- p[col + 1]

    if (end >= start) {
      rows <- i[start:end] + 1 # row indices (1-based)
      vals <- x[start:end] # values

      # append col index + value to each row's list
      for (k in seq_along(rows)) {
        r <- rows[k]
        idx_list[[r]] <- c(idx_list[[r]], col)
        val_list[[r]] <- c(val_list[[r]], vals[k])
      }
    }
  }

  # hash each row using only the sparse pattern
  out <- character(nr)
  for (r in seq_len(nr)) {
    out[r] <- digest::digest(
      list(idx_list[[r]], val_list[[r]]),
      algo = "xxhash64",
      serialize = TRUE
    )
  }

  out
}
