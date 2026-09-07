test_that("package exports the expected top-level functions", {
  expected <- c(
    "calibration.curve.fun",
    "calibration.fun",
    "compute.Q.clever.per.id",
    "estimate.alpha.fun",
    "estimate.derivative",
    "make.calibrated.contrasts",
    "prepare.initial",
    "sim.from.data",
    "sim.generic",
    "tmle.alpha.fun"
  )

  ns <- asNamespace("MultiStateTMLE")
  for (fun in expected) {
    expect_true(exists(fun, envir = ns, inherits = FALSE))
    expect_true(is.function(get(fun, envir = ns)))
  }
})
