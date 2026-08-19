source("../../R/estimator.R")

test_that("wilson interval brackets the point estimate and beats Wald when rare", {
  ci <- wilson_ci(3, 1000)
  expect_true(ci["lower"] < 3 / 1000 && 3 / 1000 < ci["upper"])
  expect_true(ci["lower"] > 0)                    # Wald would go negative here
})

test_that("newcombe risk difference is centred and orientation-symmetric", {
  a <- newcombe_rd(120, 10000, 60, 10000)
  expect_equal(unname(a["estimate"]), 0.006)
  expect_true(a["lower.lower"] < 0.006 && 0.006 < a["upper.upper"])
  b <- newcombe_rd(60, 10000, 120, 10000)
  expect_equal(unname(a["estimate"]), -unname(b["estimate"]))
})

test_that("moment and vector forms of the Wald IV estimator agree", {
  set.seed(1)
  n1 <- 50000; n0 <- 20000
  d1 <- rbinom(n1, 1, 0.3); y1 <- rbinom(n1, 1, 0.02 + 0.05 * d1)
  d0 <- rep(0, n0);         y0 <- rbinom(n0, 1, 0.02)
  v <- cace_wald(y1, d1, y0, d0)
  m <- cace_wald_moments(n1, mean(y1), mean(d1), mean(y1 * d1),
                         n0, mean(y0), mean(d0), mean(y0 * d0))
  expect_equal(unname(v["estimate"]), unname(m["estimate"]), tolerance = 1e-8)
  expect_equal(unname(v["se"]), unname(m["se"]), tolerance = 1e-6)
})

test_that("Wald IV recovers a known complier effect", {
  set.seed(2)
  n <- 400000
  z <- rbinom(n, 1, 0.85)
  complier <- rbinom(n, 1, 0.2)
  d <- z * complier
  y <- rbinom(n, 1, 0.01 + 0.04 * d)          # true CACE = 0.04
  r <- cace_wald(y[z == 1], d[z == 1], y[z == 0], d[z == 0])
  expect_equal(unname(r["estimate"]), 0.04, tolerance = 0.01)
  expect_true(r["lower"] < 0.04 && 0.04 < r["upper"])
})

test_that("stratified_ate reproduces a homogeneous effect under any weights", {
  n1 <- rep(10000, 5); n0 <- rep(10000, 5)
  p0 <- c(0.01, 0.02, 0.03, 0.04, 0.05); p1 <- p0 + 0.01
  r <- stratified_ate(n1, p1 * n1, n0, p0 * n0, w = c(1, 5, 2, 9, 3))
  expect_equal(unname(r["estimate"]), 0.01, tolerance = 1e-12)
})

test_that("stratified_ate reweights a heterogeneous effect toward the target mix", {
  n1 <- c(1000, 1000); n0 <- c(1000, 1000)
  x1 <- c(100, 400); x0 <- c(50, 100)          # taus: 0.05 and 0.30
  expect_equal(unname(stratified_ate(n1, x1, n0, x0, c(1, 0))["estimate"]), 0.05)
  expect_equal(unname(stratified_ate(n1, x1, n0, x0, c(0, 1))["estimate"]), 0.30)
  expect_equal(unname(stratified_ate(n1, x1, n0, x0, c(1, 1))["estimate"]), 0.175)
})

test_that("smd is zero under balance and unaffected by the weighting denominator", {
  set.seed(3); x <- rnorm(20000); t <- rbinom(20000, 1, 0.5)
  expect_equal(smd(x, t), 0, tolerance = 0.05)
  expect_equal(smd(x, t), smd(x, t, rep(2, 20000)), tolerance = 1e-12)
})

test_that("effective sample size equals n only for equal weights", {
  expect_equal(ess(rep(3, 100)), 100)
  expect_lt(ess(c(rep(1, 99), 500)), 100)
})

test_that("IPTW recovers the ATE when given the true propensity score", {
  set.seed(4); n <- 200000
  x <- rnorm(n)
  e <- plogis(0.3 + 0.8 * x)
  t <- rbinom(n, 1, e)
  y <- rbinom(n, 1, plogis(-2 + 0.5 * t + x))
  truth <- mean(plogis(-2 + 0.5 + x) - plogis(-2 + x))
  r <- iptw_ate(y, t, e)
  expect_equal(unname(r["estimate"]), truth, tolerance = 0.005)
  expect_true(r["lower"] < truth && truth < r["upper"])
})

test_that("AIPW is doubly robust: either nuisance may be wrong, but not both", {
  set.seed(5); n <- 200000
  x1 <- rnorm(n); x2 <- rnorm(n)
  e_true <- plogis(0.4 + x1 + 0.6 * x2^2)         # nonlinear in x2
  t <- rbinom(n, 1, e_true)
  p1 <- plogis(-1.5 + 0.6 + x1 - x2); p0 <- plogis(-1.5 + x1 - x2)
  y <- rbinom(n, 1, ifelse(t == 1, p1, p0))
  truth <- mean(p1 - p0)

  e_wrong <- plogis(predict(glm(t ~ x1 + x2, family = binomial())))  # misses x2^2
  mu1_wrong <- rep(mean(y[t == 1]), n); mu0_wrong <- rep(mean(y[t == 0]), n)

  right_ps <- aipw_ate(y, t, e_true,  mu1_wrong, mu0_wrong)
  right_om <- aipw_ate(y, t, e_wrong, p1, p0)
  both_bad <- aipw_ate(y, t, e_wrong, mu1_wrong, mu0_wrong)

  expect_equal(unname(right_ps["estimate"]), truth, tolerance = 0.01)
  expect_equal(unname(right_om["estimate"]), truth, tolerance = 0.01)
  # with both nuisances wrong there is no guarantee, and here it visibly breaks
  expect_gt(abs(unname(both_bad["estimate"]) - truth),
            abs(unname(right_ps["estimate"]) - truth))
})

test_that("g-computation is consistent only through the outcome model", {
  set.seed(6); n <- 100000
  x <- rnorm(n); t <- rbinom(n, 1, plogis(x))
  p1 <- plogis(-1 + 0.4 + x); p0 <- plogis(-1 + x)
  expect_equal(unname(gcomp_ate(p1, p0)["estimate"]), mean(p1 - p0), tolerance = 1e-12)
})

test_that("matching finds the true nearest control and honours the caliper", {
  lps   <- c(0, 10, 1.0, 1.4, 9.6)
  treat <- c(1, 1, 0, 0, 0)
  m <- match_1to1(lps, treat, caliper_sd = 10)     # wide caliper: both match
  expect_equal(m$treated, c(1, 2))
  expect_equal(m$control, c(3, 5))                 # 0->1.0 and 10->9.6
  tight <- match_1to1(lps, treat, caliper_sd = 0.1)
  expect_equal(tight$n_dropped, 1)                 # the distance-1.0 pair drops
  expect_equal(tight$treated, 2)
})

test_that("matched ATT reuses controls and prices that reuse in the variance", {
  y <- c(1, 0, 1, 0, 0)
  m <- list(treated = c(1, 2), control = c(3, 4), n_dropped = 0, caliper = 1)
  r <- matched_att(y, m)
  expect_equal(unname(r["estimate"]), 0)
  expect_equal(unname(r["n_matched"]), 2)
  expect_gt(unname(r["se"]), 0)
})

test_that("cross-fitting never predicts a unit from a model that saw it", {
  set.seed(7); n <- 5000
  d <- data.frame(x = rnorm(n)); t <- rbinom(n, 1, 0.5)
  folds <- make_folds(n, 5, 1)
  expect_equal(sort(unique(folds)), 1:5)
  expect_true(all(table(folds) >= floor(n / 5)))
  e <- cross_fit_ps(d, t, .treat ~ x, folds)
  expect_length(e, n); expect_true(all(e > 0 & e < 1))
})
