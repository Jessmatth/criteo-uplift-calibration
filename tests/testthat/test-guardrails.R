source("../../R/guardrails.R")

G <- list(max_relative_bias = 0.10, require_ci_covers_benchmark = TRUE,
          max_ci_width_ratio = 3.0, max_extreme_ps_share = 0.05,
          max_abs_smd = 0.10, min_ess_ratio = 0.10)

test_that("relative bias is scaled by the benchmark, not absolute", {
  expect_true(check_relative_bias(0.00105, 0.001, 0.10)$pass)    # 5% off
  expect_false(check_relative_bias(0.0012, 0.001, 0.10)$pass)    # 20% off
  # the same absolute error passes on a large benchmark and fails on a small one
  expect_true(check_relative_bias(0.1002, 0.100, 0.10)$pass)
  expect_false(check_relative_bias(0.0012, 0.0002, 0.10)$pass)
})

test_that("coverage is checked on the benchmark point, both directions", {
  expect_true(check_ci_covers(0.9, 1.1, 1.0, TRUE)$pass)
  expect_false(check_ci_covers(1.05, 1.2, 1.0, TRUE)$pass)
  expect_false(check_ci_covers(0.5, 0.9, 1.0, TRUE)$pass)
  expect_true(check_ci_covers(1.05, 1.2, 1.0, FALSE)$pass)       # not required
})

test_that("an estimate can be accurate and still fail on width", {
  wide <- check_ci_width(-1, 1, -0.1, 0.1, 3.0)
  expect_equal(wide$value, 10)
  expect_false(wide$pass)
  expect_true(check_ci_width(-0.2, 0.2, -0.1, 0.1, 3.0)$pass)
})

test_that("positivity, balance and ESS thresholds are one-sided as intended", {
  expect_true(check_positivity(0.02, 0.05)$pass)
  expect_false(check_positivity(0.30, 0.05)$pass)
  expect_true(check_balance(-0.08, 0.10)$pass)     # sign must not matter
  expect_false(check_balance(0.25, 0.10)$pass)
  expect_true(check_ess(0.42, 0.10)$pass)
  expect_false(check_ess(0.01, 0.10)$pass)         # weights collapsed onto few units
})

test_that("non-finite diagnostics fail rather than pass by accident", {
  expect_false(check_relative_bias(0.1, 0, 0.10)$pass)   # benchmark of zero
  expect_false(check_positivity(NaN, 0.05)$pass)
  expect_false(check_ess(NA_real_, 0.10)$pass)
})

est   <- data.frame(estimate = 0.00102, lower = 0.00090, upper = 0.00114)
bench <- data.frame(benchmark = 0.00100, benchmark_lower = 0.00094,
                    benchmark_upper = 0.00106)
good  <- data.frame(extreme_ps_share = 0.01, max_abs_smd_weighted = 0.03,
                    ess_ratio = 0.55)

test_that("a clean estimate passes every applicable check", {
  r <- evaluate_row(est, bench, good, G)
  expect_equal(r$verdict, "PASS")
  expect_equal(r$checks_applicable, 6)
  expect_equal(r$checks_failed, 0)
})

test_that("one bad diagnostic fails the row even when the estimate is close", {
  bad <- good; bad$ess_ratio <- 0.02
  r <- evaluate_row(est, bench, bad, G)
  expect_equal(r$verdict, "FAIL")
  expect_false(r$ess_pass)
  expect_true(r$relative_bias_pass)      # the point estimate was fine
})

test_that("estimators without a propensity model are scored on what applies", {
  r <- evaluate_row(est, bench, NULL, G)
  expect_equal(r$checks_applicable, 3)
  expect_true(is.na(r$positivity_pass))
  expect_equal(r$verdict, "PASS")        # NA checks must not count as failures
})

test_that("a badly biased estimate fails regardless of tidy diagnostics", {
  far <- data.frame(estimate = 0.005, lower = 0.0045, upper = 0.0055)
  r <- evaluate_row(far, bench, good, G)
  expect_equal(r$verdict, "FAIL")
  expect_false(r$relative_bias_pass)
  expect_false(r$ci_covers_pass)
})
