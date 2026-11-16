test_that("risk() and years_lost() take bare names or strings", {
  expect_equal(risk(death), risk("death"))
  expect_equal(risk(death, z)$processes, c("death", "z"))
  targets <- c("death", "z")
  expect_equal(do.call(risk, as.list(targets)), risk(death, z))
  expect_equal(years_lost(death)$estimand, "years_lost")

  expect_error(risk(), "at least one process")
  expect_error(risk(death, death), "only be given once")
  expect_error(risk(death + 1), "bare names or strings")
  expect_error(risk(c("death", "z")), "bare names or strings")
  expect_error(risk(death, foo = 1), "Unknown argument.*foo")
  expect_error(risk(death, when = z ~ z == 0), "one-sided")
  expect_error(risk(death, only_first = NA), "only_first")
})

test_that("targets and interventions print as their call", {
  expect_output(print(risk(death)), "Target: risk\\(death\\)")
  expect_equal(
    format(years_lost(death, when = ~ z == 0, only_first = TRUE)),
    "years_lost(death, when = ~z == 0, only_first = TRUE)"
  )
  expect_output(
    print(scale_intensity(z = 0.5, w = ~ time / 2)),
    "scale_intensity\\(z = 0.5, w = ~time/2\\)"
  )
  expect_equal(
    format(calibrate(risk(z), level(0.2))),
    "calibrate(risk(z), level(0.2))"
  )
  expect_output(print(relative_change(0.8)), "relative_change\\(0.8\\)")
})

test_that("scale_intensity() validates its factors", {
  expect_error(scale_intensity(), "at least one process")
  expect_error(scale_intensity(0.5), "named by process")
  expect_error(scale_intensity(z = 0.5, 2), "named by process")
  # jarl-ignore duplicated_arguments: tests the error for a duplicated process.
  expect_error(scale_intensity(z = 0.5, z = 2), "only be given once")
  expect_error(scale_intensity(z = -1), "z")
  expect_error(scale_intensity(z = "a"), "z")
  expect_error(scale_intensity(z = y ~ time), "one-sided")
})
