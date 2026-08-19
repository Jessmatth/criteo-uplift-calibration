# ---- Estimation and inference primitives ------------------------------------
# Pure functions, no I/O, unit-tested in tests/testthat/test-estimator.R.

# ---- Proportions and risk differences ---------------------------------------

wilson_ci <- function(x, n, conf = 0.95) {
  z  <- stats::qnorm(1 - (1 - conf) / 2)
  p  <- x / n
  d  <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / d
  hw  <- z / d * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  c(lower = ctr - hw, upper = ctr + hw)
}

# Newcombe's hybrid-score interval for a difference of independent proportions.
# Preferred over the Wald interval here because the conversion rate is ~0.3%:
# with a rare outcome the Wald interval is anticonservative and can cross zero
# when the score interval does not.
newcombe_rd <- function(x1, n1, x0, n0, conf = 0.95) {
  p1 <- x1 / n1; p0 <- x0 / n0
  c1 <- wilson_ci(x1, n1, conf); c0 <- wilson_ci(x0, n0, conf)
  d  <- p1 - p0
  c(estimate = d,
    lower = d - sqrt((p1 - c1["lower"])^2 + (c0["upper"] - p0)^2),
    upper = d + sqrt((c1["upper"] - p1)^2 + (p0 - c0["lower"])^2))
}

# ---- Instrumental-variable / CACE -------------------------------------------
# `exposure` is post-randomization: the ad server chose who actually saw an ad,
# using signals correlated with intent. Comparing exposed to unexposed is
# therefore not a causal contrast. Randomized assignment IS a valid instrument
# for exposure (relevance: assignment moves exposure from 0 to ~3.6%;
# exclusion: assignment can only act through being shown the ad; monotonicity:
# no user is exposed only when withheld, which holds by construction here since
# exposure is impossible under control). The Wald ratio then identifies the
# effect on compliers - the users whose exposure the campaign actually changed.
cace_wald <- function(y1, d1, y0, d0) {
  n1 <- length(y1); n0 <- length(y0)
  a1 <- mean(y1); a0 <- mean(y0); b1 <- mean(d1); b0 <- mean(d0)
  num <- a1 - a0; den <- b1 - b0
  est <- num / den

  # Delta method. Y and D are correlated within an arm, so the covariance term
  # is retained; dropping it materially understates the standard error.
  S1 <- stats::cov(cbind(y1, d1)) / n1
  S0 <- stats::cov(cbind(y0, d0)) / n0
  g  <- c(1 / den, -num / den^2)                       # d(est)/d(a1,b1)
  h  <- c(-1 / den,  num / den^2)                      # d(est)/d(a0,b0)
  v  <- as.numeric(t(g) %*% S1 %*% g + t(h) %*% S0 %*% h)
  se <- sqrt(v)
  c(estimate = est, se = se,
    lower = est - 1.96 * se, upper = est + 1.96 * se,
    compliance = den)
}

# Moment form of the same estimator, so the benchmark can be computed from SQL
# aggregates over 14M rows instead of pulling them into memory. y and d are
# binary, so var = p(1-p) and cov = E[yd] - E[y]E[d].
cace_wald_moments <- function(n1, ybar1, dbar1, ydbar1, n0, ybar0, dbar0, ydbar0) {
  num <- ybar1 - ybar0; den <- dbar1 - dbar0
  est <- num / den
  S1 <- matrix(c(ybar1 * (1 - ybar1), ydbar1 - ybar1 * dbar1,
                 ydbar1 - ybar1 * dbar1, dbar1 * (1 - dbar1)), 2, 2) / n1
  S0 <- matrix(c(ybar0 * (1 - ybar0), ydbar0 - ybar0 * dbar0,
                 ydbar0 - ybar0 * dbar0, dbar0 * (1 - dbar0)), 2, 2) / n0
  g <- c(1 / den, -num / den^2); h <- c(-1 / den, num / den^2)
  se <- sqrt(as.numeric(t(g) %*% S1 %*% g + t(h) %*% S0 %*% h))
  c(estimate = est, se = se, lower = est - 1.96 * se, upper = est + 1.96 * se,
    compliance = den)
}

# Baseline conversion rate of compliers, backed out of the mixture in the
# control arm. Compliers are the users the ad server chose to reach; if their
# untreated conversion rate is far above the never-takers', that gap IS the
# targeting selection, measured rather than assumed.
complier_baseline <- function(p_control, p_treated_unexposed, compliance) {
  (p_control - (1 - compliance) * p_treated_unexposed) / compliance
}

# ---- Stratified randomized ground truth -------------------------------------
# Within a stratum, treatment is still randomized, so the difference in means is
# unbiased with no modelling. Retention into the observational sample is
# constant within a stratum, which leaves the within-stratum covariate
# distribution untouched; the benchmark for that sample is therefore the
# stratum-specific randomized effects reweighted to the sample's stratum shares.
stratified_ate <- function(n1, x1, n0, x0, w) {
  keep <- n1 > 0 & n0 > 0 & w > 0
  n1 <- n1[keep]; x1 <- x1[keep]; n0 <- n0[keep]; x0 <- x0[keep]; w <- w[keep]
  w <- w / sum(w)
  p1 <- x1 / n1; p0 <- x0 / n0
  tau_k <- p1 - p0
  var_k <- p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0
  est <- sum(w * tau_k)
  se  <- sqrt(sum(w^2 * var_k))
  c(estimate = est, se = se, lower = est - 1.96 * se, upper = est + 1.96 * se)
}

# ---- Balance diagnostics ----------------------------------------------------

weighted_mean <- function(x, w) sum(w * x) / sum(w)
weighted_var  <- function(x, w) {
  m  <- weighted_mean(x, w)
  sw <- sum(w)
  sum(w * (x - m)^2) / (sw - sum(w^2) / sw)
}

# Standardized mean difference. The denominator is pooled over the UNWEIGHTED
# groups and held fixed across adjustment methods, so that a method cannot
# appear to improve balance merely by inflating the variance it divides by.
smd <- function(x, treat, w = NULL) {
  if (is.null(w)) w <- rep(1, length(x))
  t1 <- treat == 1
  denom <- sqrt((stats::var(x[t1]) + stats::var(x[!t1])) / 2)
  if (!is.finite(denom) || denom == 0) return(0)
  (weighted_mean(x[t1], w[t1]) - weighted_mean(x[!t1], w[!t1])) / denom
}

# Kish effective sample size: how many independent observations the weighted
# sample is actually worth. A weighting scheme that balances covariates by
# concentrating on a handful of units is not a usable measurement instrument.
ess <- function(w) sum(w)^2 / sum(w^2)

# ---- Propensity scores ------------------------------------------------------

truncate_ps <- function(e, bounds) pmin(pmax(e, bounds[1]), bounds[2])

ate_weights <- function(treat, e) ifelse(treat == 1, 1 / e, 1 / (1 - e))

# Stabilized ATE weights: numerator is the marginal treatment probability.
# Stabilization does not change the estimand, it shrinks the weight variance,
# which is what keeps the CI width guardrail attainable.
stabilized_weights <- function(treat, e) {
  pt <- mean(treat)
  ifelse(treat == 1, pt / e, (1 - pt) / (1 - e))
}

# ---- IPTW -------------------------------------------------------------------
# Hajek (self-normalized) form. The influence-function SE treats the propensity
# score as known, which is conservative: estimating it would remove variance.
iptw_ate <- function(y, treat, e) {
  w  <- stabilized_weights(treat, e)
  t1 <- treat == 1
  mu1 <- weighted_mean(y[t1], w[t1]); mu0 <- weighted_mean(y[!t1], w[!t1])
  est <- mu1 - mu0
  n <- length(y)
  psi <- (treat / e) * (y - mu1) / mean(treat / e) -
         ((1 - treat) / (1 - e)) * (y - mu0) / mean((1 - treat) / (1 - e))
  se <- stats::sd(psi) / sqrt(n)
  c(estimate = est, se = se, lower = est - 1.96 * se, upper = est + 1.96 * se,
    ess_ratio = ess(w) / n)
}

# ---- AIPW / doubly robust ---------------------------------------------------
# Consistent if EITHER the propensity model or the outcome model is right. The
# nuisances must be fit out-of-fold (see cross_fit_nuisance) or the influence
# function SE understates the true variance through overfitting.
aipw_ate <- function(y, treat, e, mu1, mu0) {
  psi <- mu1 - mu0 +
    treat * (y - mu1) / e -
    (1 - treat) * (y - mu0) / (1 - e)
  est <- mean(psi)
  se  <- stats::sd(psi) / sqrt(length(psi))
  c(estimate = est, se = se, lower = est - 1.96 * se, upper = est + 1.96 * se)
}

# ---- Matching ---------------------------------------------------------------
# 1:1 nearest neighbour on the LINEAR propensity score, with replacement, inside
# a caliper. Matching on the linear score rather than the probability keeps
# distances meaningful in the tails, where a rare outcome puts most of the mass.
#
# With replacement is chosen deliberately: without replacement the estimate
# depends on the order units are processed, and with 11.6% of rows sharing an
# identical feature vector that ordering is close to arbitrary.
#
# Implementation is a sorted binary search rather than a pairwise distance
# matrix: at 2M rows the latter is not computable.
match_1to1 <- function(lps, treat, caliper_sd = 0.2) {
  cal <- caliper_sd * stats::sd(lps)
  it  <- which(treat == 1); ic <- which(treat == 0)
  cs  <- lps[ic]; ord <- order(cs); cs_sorted <- cs[ord]; ic_sorted <- ic[ord]

  pos <- findInterval(lps[it], cs_sorted)
  lo  <- pmax(pos, 1L); hi <- pmin(pos + 1L, length(cs_sorted))
  d_lo <- abs(lps[it] - cs_sorted[lo]); d_hi <- abs(lps[it] - cs_sorted[hi])
  pick <- ifelse(d_lo <= d_hi, lo, hi)
  dist <- pmin(d_lo, d_hi)

  ok <- dist <= cal
  list(treated = it[ok], control = ic_sorted[pick[ok]],
       n_dropped = sum(!ok), caliper = cal)
}

# ATT from a matched set. Variance conditions on the matched sets and accounts
# for control reuse through K_i, the number of times a control is drawn.
# Abadie & Imbens (2008) show the bootstrap is INVALID for nearest-neighbour
# matching estimators, so it is not offered here.
matched_att <- function(y, m) {
  yt <- y[m$treated]; yc <- y[m$control]
  n1 <- length(yt)
  est <- mean(yt) - mean(yc)
  K  <- table(m$control)
  yc_u <- y[as.integer(names(K))]
  v <- (sum((yt - mean(yt))^2) + sum(as.numeric(K)^2 * (yc_u - mean(yc))^2)) / n1^2
  se <- sqrt(v)
  c(estimate = est, se = se, lower = est - 1.96 * se, upper = est + 1.96 * se,
    n_matched = n1, n_dropped = m$n_dropped)
}

# ---- Cross-fitting ----------------------------------------------------------
# Sample splitting for AIPW. Fitting nuisances on the same rows used for the
# estimate reintroduces the overfitting bias that double robustness is supposed
# to remove; K-fold cross-fitting restores the root-n behaviour.
make_folds <- function(n, k, seed) {
  set.seed(seed)
  sample(rep_len(seq_len(k), n))
}

# Nuisances are cross-fit separately so that a propensity spec and an outcome
# spec can be paired without refitting either. Every prediction a unit receives
# comes from a model that never saw that unit.
cross_fit_ps <- function(data, treat, formula, folds) {
  e <- numeric(nrow(data))
  for (k in sort(unique(folds))) {
    tr <- folds != k; te <- folds == k
    fit <- stats::glm(formula, family = stats::binomial(),
                      data = cbind(data[tr, , drop = FALSE], .treat = treat[tr]))
    e[te] <- stats::predict(fit, newdata = data[te, , drop = FALSE], type = "response")
  }
  e
}

cross_fit_om <- function(data, y, treat, formula, folds) {
  mu1 <- numeric(nrow(data)); mu0 <- numeric(nrow(data))
  for (k in sort(unique(folds))) {
    tr <- folds != k; te <- folds == k
    for (a in c(0, 1)) {
      idx <- which(tr & treat == a)
      fit <- stats::glm(formula, family = stats::binomial(),
                        data = cbind(data[idx, , drop = FALSE], .y = y[idx]))
      p <- stats::predict(fit, newdata = data[te, , drop = FALSE], type = "response")
      if (a == 1) mu1[te] <- p else mu0[te] <- p
    }
  }
  list(mu1 = mu1, mu0 = mu0)
}

# Plug-in g-computation: average the two modelled outcome surfaces over the
# sample. Consistent only if the outcome model is right - it has no second
# chance, which is exactly what makes it the useful contrast for AIPW.
gcomp_ate <- function(mu1, mu0) {
  psi <- mu1 - mu0
  est <- mean(psi); se <- stats::sd(psi) / sqrt(length(psi))
  c(estimate = est, se = se, lower = est - 1.96 * se, upper = est + 1.96 * se)
}
