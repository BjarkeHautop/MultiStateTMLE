### script.t2d.R ---
#----------------------------------------------------------------------
## Author: Helene
## Created: Sep  3 2026 (07:41)
## Version:
## Last-Updated: Sep  4 2026 (10:31)
##           By: Helene
##     Update #: 176
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:

library(hdnom)
library(MASS)
library(SuperLearner)
library(data.table)
library(survival)
library(prodlim)
library(zoo)
library(nleqslv)
library(foreach)
library(doParallel)
library(parallel)
library(xtable)
library(stringr)
library(glmnet)
library(Matrix)
library(digest)
library(rlist)
library(devtools)
#remotes::install_github("miclukacova/simevent")
library(simevent)

#-------------------------------------------------------------------------------------------#
## source relevant scripts from recurrent repository
#-------------------------------------------------------------------------------------------#

source("./usage/master.R")
run.master()

summary <- readRDS(file = paste0(".output/", "summary_for_helene", ".rds"))

set.seed(15)

dsim <- sim.from.data(n = summary$n.subjects, summary, browse = FALSE)

dsim[, summary(time)]
dsim[, table(delta)]

dsim[, table(t2d)]
dsim[, table(hf)]

process.names <- summary$model.structure$process.names
process.deltas <- summary$model.structure$process.deltas

process.names.to.deltas <- seq_along(process.names[order(process.deltas)]) - 1L
names(process.names.to.deltas) <- process.names[order(process.deltas)]

baseline.model <- paste0(names(summary$baseline.summary), collapse = "+")

tau <- 1250

Nsim <- 2500
use.cores <- 1

dsim <- sim.from.data(n = Nsim, summary, browse = FALSE)

## fit initial, here based on same models as Jens

init.fit <- prepare.initial(
  dsim,
  use.cores = use.cores,
  tau = tau,
  fit.types = list(
    outcome1 = list(
      model = paste0(
        "Surv(tstart, tstop, delta == ",
        process.names.to.deltas["outcome1"],
        ")~t2d+hf+",
        baseline.model
      ),
      fit = "cox"
    ),

    cr1 = list(
      model = paste0(
        "Surv(tstart, tstop, delta == ",
        process.names.to.deltas["cr1"],
        ")~t2d+hf+",
        baseline.model
      ),
      fit = "cox"
    ),

    censoring = list(
      model = paste0(
        "Surv(tstart, tstop, delta == ",
        process.names.to.deltas["censoring"],
        ")~t2d+hf+",
        baseline.model
      ),
      fit = "cox"
    ),

    t2d = list(
      model = paste0(
        "Surv(tstart, tstop, delta == ",
        process.names.to.deltas["t2d"],
        ")~hf+",
        baseline.model
      ),
      fit = "cox",
      at.risk = function(dt) dt[["t2d"]] == 0
    ),

    hf = list(
      model = paste0(
        "Surv(tstart, tstop, delta == ",
        process.names.to.deltas["hf"],
        ")~t2d+",
        baseline.model
      ),
      fit = "cox",
      at.risk = function(dt) dt[["hf"]] == 0
    )
  ),
  fit.treatment = paste0(
    "arm_n~",
    gsub("\\+arm_n|arm_n\\+", "", baseline.model)
  ),
  a = 0:1,
  verbose = TRUE
)

###################################################
## estimation of calibrated target parameters

if (FALSE) {
  ## (((we don't use this, but we could have)))
  est.placebo.at.sema.burden <- calibration.fun(
    initial.fit = init.fit,
    target = c("outcome1", "cr1"),
    z.name = "t2d",
    a = 0,
    rho.1a = 1,
    tau = tau,
    use.cores = use.cores,
    output.eic = TRUE
  )
}

est.sema.at.placebo.burden <- calibration.fun(
  initial.fit = init.fit,
  target = c("outcome1", "cr1"),
  z.name = "t2d",
  a = 1,
  rho.1a = 1,
  tau = tau,
  use.cores = use.cores,
  output.eic = TRUE
)

est.outcome1.sema <- tmle.alpha.fun(
  initial.fit = init.fit,
  target = "outcome1",
  z.name = "t2d",
  tau = tau,
  alpha = 1,
  a = 1,
  use.cores = use.cores,
  output.weights = TRUE,
  output.eic = TRUE,
  verbose = TRUE
)


est.cr1.sema <- tmle.alpha.fun(
  initial.fit = init.fit,
  target = "cr1",
  z.name = "t2d",
  tau = tau,
  alpha = 1,
  a = 1,
  use.cores = use.cores,
  output.weights = TRUE,
  output.eic = TRUE,
  verbose = TRUE
)

est.outcome1.placebo <- tmle.alpha.fun(
  initial.fit = init.fit,
  target = "outcome1",
  z.name = "t2d",
  tau = tau,
  alpha = 1,
  a = 0,
  use.cores = use.cores,
  output.weights = TRUE,
  output.eic = TRUE,
  verbose = TRUE
)

est.cr1.placebo <- tmle.alpha.fun(
  initial.fit = init.fit,
  target = "cr1",
  z.name = "t2d",
  tau = tau,
  alpha = 1,
  a = 0,
  use.cores = use.cores,
  output.weights = TRUE,
  output.eic = TRUE,
  verbose = TRUE
)

(outcome1.calibration.contrasts <-
  make.calibrated.contrasts(
    target = "outcome1",
    treatment.fit = est.outcome1.sema,
    placebo.fit = est.outcome1.placebo,
    calibrated.fit = est.sema.at.placebo.burden
  ))

(cr1.calibration.contrasts <-
  make.calibrated.contrasts(
    target = "cr1",
    treatment.fit = est.cr1.sema,
    placebo.fit = est.cr1.placebo,
    calibrated.fit = est.sema.at.placebo.burden
  ))

if (FALSE) {
  ### doing the same as above, but manually

  est.outcome1.sema$estimate[["tmle.est"]]
  est.sema.at.placebo.burden$estimate[["target.est.outcome1"]]

  est.cr1.sema$estimate[["tmle.est"]]
  est.sema.at.placebo.burden$estimate[["target.est.cr1"]]

  (outcome1.total <-
    est.outcome1.sema$estimate[["tmle.est"]] -
    est.outcome1.placebo$estimate[["tmle.est"]])

  (outcome1.total.se <-
    sqrt(mean(
      (est.outcome1.sema$eic -
        est.outcome1.placebo$eic)^2 /
        length(est.outcome1.sema$eic)
    )))

  (outcome1.mediated <-
    est.outcome1.sema$estimate[["tmle.est"]] -
    est.sema.at.placebo.burden$estimate[["target.est.outcome1"]])

  (outcome1.mediated.se <-
    sqrt(mean(
      (est.outcome1.sema$eic -
        est.sema.at.placebo.burden$eic$outcome1)^2 /
        length(est.outcome1.sema$eic)
    )))

  (outcome1.calibrated <-
    est.sema.at.placebo.burden$estimate[["target.est.outcome1"]] -
    est.outcome1.placebo$estimate[["tmle.est"]])

  (outcome1.calibrated.se <-
    sqrt(mean(
      (est.sema.at.placebo.burden$eic$outcome1 -
        est.outcome1.placebo$eic)^2 /
        length(est.outcome1.sema$eic)
    )))

  (cr1.total <-
    est.cr1.sema$estimate[["tmle.est"]] -
    est.cr1.placebo$estimate[["tmle.est"]])

  (cr1.total.se <-
    sqrt(mean(
      (est.cr1.sema$eic -
        est.cr1.placebo$eic)^2 /
        length(est.cr1.sema$eic)
    )))

  (cr1.mediated <-
    est.cr1.sema$estimate[["tmle.est"]] -
    est.sema.at.placebo.burden$estimate[["target.est.cr1"]])

  (cr1.mediated.se <-
    sqrt(mean(
      (est.cr1.sema$eic -
        est.sema.at.placebo.burden$eic$cr1)^2 /
        length(est.cr1.sema$eic)
    )))

  (cr1.calibrated <-
    est.sema.at.placebo.burden$estimate[["target.est.cr1"]] -
    est.cr1.placebo$estimate[["tmle.est"]])

  (cr1.calibrated.se <-
    sqrt(mean(
      (est.sema.at.placebo.burden$eic$cr1 -
        est.cr1.placebo$eic)^2 /
        length(est.cr1.sema$eic)
    )))

  (outcome1.estimates <- c(
    total = outcome1.total,
    mediated = outcome1.mediated,
    calibrated = outcome1.calibrated,
    decomposition.error = outcome1.total -
      outcome1.calibrated -
      outcome1.mediated
  ))

  outcome1.total + c(-1, 1) * 1.96 * outcome1.total.se
  outcome1.mediated + c(-1, 1) * 1.96 * outcome1.mediated.se
  outcome1.calibrated + c(-1, 1) * 1.96 * outcome1.calibrated.se

  (cr1.estimates <- c(
    total = cr1.total,
    mediated = cr1.mediated,
    calibrated = cr1.calibrated,
    decomposition.error = cr1.total -
      cr1.mediated -
      cr1.calibrated
  ))

  cr1.total + c(-1, 1) * 1.96 * cr1.total.se
  cr1.mediated + c(-1, 1) * 1.96 * cr1.mediated.se
  cr1.calibrated + c(-1, 1) * 1.96 * cr1.calibrated.se
}

###################################################
## estimation of target parameters, varying alphas

if (FALSE) {
  use.cores <- 30

  est.calibration.curve.a1 <- calibration.curve.fun(
    initial.fit = init.fit,
    target = c("outcome1", "cr1"),
    z.name = "t2d",
    alpha.grid = seq(0, 10, length = 50),
    a = 1,
    tau = tau,
    use.cores = use.cores,
    output.eic = TRUE
  )

  saveRDS(
    est.calibration.curve.a1,
    file = paste0(
      "./output/",
      "simulated-t2d-tau1250-",
      "est-calibration-curve-a1",
      ".rds"
    )
  )

  est.calibration.curve.a0 <- calibration.curve.fun(
    initial.fit = init.fit,
    target = c("outcome1", "cr1"),
    z.name = "t2d",
    alpha.grid = seq(0, 10, length = 50),
    a = 0,
    tau = tau,
    use.cores = use.cores,
    output.eic = TRUE
  )

  saveRDS(
    est.calibration.curve.a0,
    file = paste0(
      "./output/",
      "simulated-t2d-tau1250-",
      "est-calibration-curve-a0",
      ".rds"
    )
  )
}

###################################################
## estimation of target parameters, varying alphas

if (FALSE) {
  tmle.est.list <- list()

  for (aa in 0:1) {
    print(aa)
    for (alpha in c(0.25, 1, 1.5)) {
      print(alpha)
      tmle.est.list[[paste0("treat = ", aa)]][[paste0("alpha = ", alpha)]][[
        "outcome1"
      ]] <-
        tmle.alpha.fun(
          initial.fit = init.fit,
          target = "outcome1",
          z.name = "t2d",
          tau = tau,
          alpha = alpha,
          a = aa,
          target.by.state = TRUE,
          use.cores = use.cores,
          output.weights = TRUE,
          output.eic = FALSE,
          verbose = TRUE
        )
      tmle.est.list[[paste0("treat = ", aa)]][[paste0("alpha = ", alpha)]][[
        "cr1"
      ]] <-
        tmle.alpha.fun(
          initial.fit = init.fit,
          target = "cr1",
          z.name = "t2d",
          tau = tau,
          alpha = alpha,
          a = aa,
          target.by.state = TRUE,
          use.cores = use.cores,
          output.weights = TRUE,
          output.eic = FALSE,
          verbose = TRUE
        )
      tmle.est.list[[paste0("m=", m)]][[paste0("treat = ", aa)]][[paste0(
        "alpha = ",
        alpha
      )]][["t2d"]] <-
        tmle.alpha.fun(
          initial.fit = init.fit,
          target = "t2d",
          z.name = "t2d",
          tau = tau,
          alpha = alpha,
          a = aa,
          target.by.state = TRUE,
          use.cores = use.cores,
          output.weights = TRUE,
          output.eic = FALSE,
          verbose = TRUE
        )
    }
  }
}

######################################################################
### script.t2d.R ends here
