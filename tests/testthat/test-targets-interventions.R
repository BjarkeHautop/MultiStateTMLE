test_that("risk() and years_lost() take bare names or strings", {
  expect_equal(risk(death), risk("death"))
  expect_equal(risk(death, illness)$processes, c("death", "illness"))
  targets <- c("death", "illness")
  expect_equal(do.call(risk, as.list(targets)), risk(death, illness))
  expect_equal(years_lost(death)$estimand, "years_lost")

  expect_error(risk(), "at least one process")
  expect_error(risk(death, death), "only be given once")
  expect_error(risk(death + 1), "bare names or strings")
  expect_error(risk(c("death", "illness")), "bare names or strings")
  expect_error(risk(death, foo = 1), "Unknown argument.*foo")
  expect_error(risk(death, when = illness ~ illness == 0), "one-sided")
  expect_error(risk(death, only_first = NA), "only_first")
})

test_that("targets and interventions print as their call", {
  expect_output(print(risk(death)), "Target: risk\\(death\\)")
  expect_equal(
    format(years_lost(death, when = ~ illness == 0, only_first = TRUE)),
    "years_lost(death, when = ~illness == 0, only_first = TRUE)"
  )
  expect_output(
    print(scale_intensity(illness = 0.5, w = ~ time / 2)),
    "scale_intensity\\(illness = 0.5, w = ~time/2\\)"
  )
  expect_equal(
    format(calibrate(risk(illness), level(0.2))),
    "calibrate(risk(illness), level(0.2))"
  )
  expect_output(print(relative_change(0.8)), "relative_change\\(0.8\\)")
})

test_that("scale_intensity() validates its factors", {
  expect_error(scale_intensity(), "at least one process")
  expect_error(scale_intensity(0.5), "named by process")
  expect_error(scale_intensity(illness = 0.5, 2), "named by process")
  # jarl-ignore duplicated_arguments: tests the error for a duplicated process.
  expect_error(
    scale_intensity(illness = 0.5, illness = 2),
    "only be given once"
  )
  expect_error(scale_intensity(illness = -1), "illness")
  expect_error(scale_intensity(illness = "a"), "illness")
  expect_error(scale_intensity(illness = y ~ time), "one-sided")
})
