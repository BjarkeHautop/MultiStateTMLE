# Simulate the `illness_death` example dataset with simevent: a trial
# where a randomised treatment A0 lowers the hazard of illness (at most
# once), which raises the hazard of death. Treatment also lowers the death
# hazard directly, so its effect on death is partly through the illness. The
# illness hazard depends on a baseline covariate L0, and follow-up is subject
# to independent censoring.

library(simevent)

set.seed(1405)
model <- sim_model(
  L0 = sim_covariate(function(N) rnorm(N)),
  A0 = sim_covariate(function(N) rbinom(N, 1, 0.5)),
  illness = sim_process("transient", limit = 1, eta = 0.7, nu = 1),
  death = sim_process("terminal", eta = 0.3, nu = 1),
  censoring = sim_process("censoring", eta = 0.1, nu = 1),
  effects = list(
    sim_effect(illness ~ L0 + A0, coef = c(0.5, -1)),
    sim_effect(death ~ illness + A0, coef = c(1.2, -0.3))
  )
)
sim <- simulate_events(model, n = 200)

illness_death <- data.frame(
  id = sim$id,
  time = sim$time,
  event = factor(sim$event, levels = c("illness", "death", "censoring")),
  L0 = sim$L0,
  A0 = as.integer(sim$A0)
)

usethis::use_data(illness_death, overwrite = TRUE)
