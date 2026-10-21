# Simulate the `illness_death` example dataset with simevent: an illness z
# (at most once) that raises the hazard of death death, both depending on
# a baseline covariate L0, with independent censoring.

library(simevent)
library(data.table)

set.seed(1405)
model <- sim_model(
  L0 = sim_covariate(function(N) rnorm(N)),
  z = sim_process("transient", limit = 1, eta = 0.2, nu = 1),
  death = sim_process("terminal", eta = 0.3, nu = 1),
  censoring = sim_process("censoring", eta = 0.1, nu = 1),
  effects = list(
    sim_effect("L0", "z", 0.5),
    sim_effect("z", "death", 0.7)
  )
)
sim <- sim_events(model, n = 150)

codes <- c(censoring = 0L, death = 1L, z = 2L)
illness_death <- data.frame(
  id = sim$id,
  time = sim$time,
  delta = unname(codes[as.character(sim$event)]),
  L0 = sim$L0
)

usethis::use_data(illness_death, overwrite = TRUE)
