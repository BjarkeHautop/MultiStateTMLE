# Shared helpers for tests of recurrent processes: simulate from a known
# DGP where relapse is recurrent and the death hazard grows with the number
# of relapses so far, and compute the truth in closed form.
#
# DGP, all hazards constant in time:
# * relapse: 0.5 * exp(0.5 * biomarker), any number of times,
# * death: 0.2 * exp(0.4 * k) after k relapses, ends follow-up,
# * censoring: 0.1, ends follow-up,
# with biomarker ~ N(0, 1).

simulate.recurrent <- function(n, seed) {
  model <- simevent::sim_model(
    biomarker = simevent::sim_covariate(function(N) rnorm(N)),
    relapse = simevent::sim_process("transient", eta = 0.5, nu = 1),
    death = simevent::sim_process("terminal", eta = 0.2, nu = 1),
    censoring = simevent::sim_process("censoring", eta = 0.1, nu = 1),
    effects = list(
      simevent::sim_effect(relapse ~ biomarker, coef = 0.5),
      simevent::sim_effect(death ~ relapse, coef = 0.4)
    )
  )
  simevent::simulate_events(model, n = n, seed = seed)
}

recurrent.models <- list(
  relapse ~ biomarker,
  death ~ biomarker + relapse,
  censoring ~ biomarker
)

recurrent.fit <- function(n = 80, seed = 1, max_count = 3) {
  fit_initial(
    simulate.recurrent(n, seed),
    models = recurrent.models,
    tau = 1,
    max_count = max_count
  )
}

# The true targets by tau = 1 when the relapse hazard after k relapses, in
# the time interval number `piece` of `breaks`, is scaled by
# alpha(k, piece). Given biomarker, the process is a Markov chain on (alive or
# dead) x (relapses so far, truncated at K), whose transition probabilities
# and expected occupation times solve the forward equations; these are
# integrated over biomarker ~ N(0, 1) numerically.
recurrent.truth <- function(alpha = function(k, piece) 1, breaks = 1) {
  by_L0 <- function(biomarker, K = 15) {
    nk <- K + 1
    S <- 2 * nk
    alive <- 1:nk
    dead <- nk + 1:nk
    h_relapse <- 0.5 * exp(0.5 * biomarker)
    p <- c(1, rep(0, S - 1))
    occupation <- rep(0, S)
    n_relapse <- 0
    start <- 0
    for (piece in seq_along(breaks)) {
      Q <- matrix(0, S, S)
      rate <- numeric(nk)
      for (k in 0:K) {
        i <- k + 1
        rate[i] <- if (k < K) alpha(k, piece) * h_relapse else 0
        if (k < K) {
          Q[i, i + 1] <- rate[i]
        }
        Q[i, nk + i] <- 0.2 * exp(0.4 * k)
        Q[i, i] <- -sum(Q[i, -i])
      }
      # The top-right block of exp(t * [[Q, I], [0, 0]]) holds the
      # integrals of exp(s * Q) over s in [0, t], the occupation times.
      A <- rbind(cbind(Q, diag(S)), matrix(0, S, 2 * S))
      E <- as.matrix(Matrix::expm(Matrix::Matrix(A * (breaks[piece] - start))))
      occupied <- drop(p %*% E[1:S, S + 1:S])
      n_relapse <- n_relapse + sum(rate * occupied[alive])
      occupation <- occupation + occupied
      p <- drop(p %*% E[1:S, 1:S])
      start <- breaks[piece]
    }
    c(
      death = sum(p[dead]),
      death_before_relapse = p[dead[1]],
      relapse = n_relapse,
      relapse_first = 1 - p[alive[1]] - p[dead[1]],
      years_lost_death = sum(occupation[dead])
    )
  }
  nodes <- seq(-8, 8, length.out = 401)
  colSums(t(sapply(nodes, by_L0)) * dnorm(nodes) * diff(nodes[1:2]))
}

# The truth under a constant scaling of every relapse.
recurrent.truth.constant <- function(alpha) {
  recurrent.truth(function(k, piece) alpha)
}
