source("../../R/confound.R")

test_that("solve_intercept hits the requested sample size at every alpha", {
  n_k <- rep(100000, 20); m_k <- scale(seq_len(20))[, 1]
  for (alpha in c(0, 1, 2, 4)) {
    a <- solve_intercept(n_k, m_k, alpha, +1, target_n = 300000)
    expect_equal(sum(n_k * retention(a, m_k, alpha, +1)), 300000, tolerance = 1e-4)
  }
})

test_that("holding n fixed across alpha is what makes scenarios comparable", {
  n_k <- rep(50000, 10); m_k <- scale(seq_len(10))[, 1]
  sizes <- sapply(c(0, 1, 2, 4), function(al)
    sum(n_k * retention(solve_intercept(n_k, m_k, al, +1, 1e5), m_k, al, +1)))
  expect_equal(sizes, rep(1e5, 4), tolerance = 1e-4)
})

test_that("alpha=0 induces no selection and larger alpha induces more", {
  n_k <- rep(50000, 10); m_k <- scale(seq_len(10))[, 1]
  r0 <- retention(solve_intercept(n_k, m_k, 0, +1, 1e5), m_k, 0, +1)
  expect_equal(diff(range(r0)), 0, tolerance = 1e-9)
  ratio <- sapply(c(1, 2, 4), function(al) {
    r <- retention(solve_intercept(n_k, m_k, al, +1, 1e5), m_k, al, +1)
    max(r) / min(r)
  })
  expect_true(all(diff(ratio) > 0))
})

test_that("treated and control selection push in opposite directions", {
  n_k <- rep(50000, 10); m_k <- scale(seq_len(10))[, 1]
  rt <- retention(solve_intercept(n_k, m_k, 2, +1, 1e5), m_k, 2, +1)
  rc <- retention(solve_intercept(n_k, m_k, 2, -1, 1e5), m_k, 2, -1)
  expect_gt(rt[10], rt[1])     # treated concentrate in high-propensity strata
  expect_lt(rc[10], rc[1])     # controls are pushed out of them
})

test_that("an unreachable target size fails loudly rather than silently clipping", {
  expect_error(solve_intercept(rep(1000, 5), scale(1:5)[, 1], 1, +1, 999999),
               "unreachable")
})

test_that("scenario ids are stable and the sample query is scenario-specific", {
  expect_equal(scenario_id("hidden", 2), "hidden_a2")
  q <- obs_sample_sql("observed", 1, 42, 50)
  expect_true(grepl("stratum_observed", q))
  expect_true(grepl("s.arm = 'observed'", q, fixed = TRUE))
  expect_true(grepl("p_treated", q))
})
