test_that("mirai daemons give the same results as running sequentially", {
  skip_on_cran()
  # The daemons run the installed package, so only test it when that is the
  # package being tested (as in R CMD check), not with devtools::load_all().
  installed <- dir.exists(
    file.path(getNamespaceInfo("MultiStateTMLE", "path"), "Meta")
  )
  skip_if_not(installed, "MultiStateTMLE is not loaded from an installation")

  models <- list(
    hazard_model(illness ~ biomarker + treatment, limit = 1),
    death ~ biomarker + treatment + illness,
    censoring ~ biomarker
  )
  run <- function() {
    fit <- fit_initial(
      fixture.data(60),
      models,
      tau = 1,
      treatment_model = treatment ~ biomarker
    )
    estimate <- tmle_estimate(
      fit,
      risk(death),
      scale_intensity(illness = 0.5),
      a = 1
    )
    list(fit = fit, estimate = estimate)
  }
  sequential <- run()

  mirai::daemons(2)
  withr::defer(mirai::daemons(0))
  expect_equal(.n_workers(), 2L)
  parallel <- run()

  expect_equal(parallel$fit$tmp_long, sequential$fit$tmp_long)
  expect_equal(parallel$fit$depend_matrix, sequential$fit$depend_matrix)
  expect_equal(
    parallel$estimate[c("estimate", "eic")],
    sequential$estimate[c("estimate", "eic")]
  )
})

test_that(".chunks() splits into runs of consecutive elements", {
  expect_equal(.chunks(1:5, 2), list(1:2, 3:5))
  expect_equal(.chunks(1:3, 8), list(1L, 2L, 3L))
  expect_equal(.chunks(7L, 4), list(7L))
  expect_equal(unlist(.chunks(1:100, 7)), 1:100)
})
