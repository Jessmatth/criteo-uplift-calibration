# ---- Confounding construction helpers ---------------------------------------
# Pure functions, unit-tested in tests/testthat/test-confound.R.

expit <- function(x) 1 / (1 + exp(-x))

# Choose the intercept that makes stratum-varying retention hit an exact target
# count. Retention is expit(a + sign * alpha * m_k); alpha sets how steeply the
# ad server favours high-propensity strata, and `a` is solved so that the sample
# size does not change as alpha does. Holding n fixed across the sweep is what
# makes bias comparable across scenarios: any movement is confounding, not
# precision.
solve_intercept <- function(n_k, m_k, alpha, sign, target_n) {
  f <- function(a) sum(n_k * expit(a + sign * alpha * m_k)) - target_n
  lo <- -50; hi <- 50
  if (f(lo) > 0 || f(hi) < 0)
    stop("target sample size unreachable at alpha=", alpha,
         "; reduce N_OBS or the alpha grid")
  stats::uniroot(f, c(lo, hi), tol = 1e-10)$root
}

retention <- function(a, m_k, alpha, sign) expit(a + sign * alpha * m_k)

scenario_id <- function(arm, alpha) sprintf("%s_a%s", arm, format(alpha))

# Deterministic pseudo-random retention: a user's draw is a hash of their id and
# the scenario seed, so a scenario is reproducible without storing 2M row ids
# per cell, and re-running does not quietly resample.
obs_sample_sql <- function(arm, alpha, seed, k_strata) {
  stratum <- sprintf("stratum_%s", arm)
  sprintf("
    WITH draw AS (
      SELECT c.*, t.%1$s AS stratum,
             (hash(c.user_id * 1000003 + %2$d) %% 1000000) / 1000000.0 AS u
      FROM criteo c JOIN targeting t USING (user_id)
    )
    SELECT d.* EXCLUDE (u), s.e_true
    FROM draw d
    JOIN scenarios s
      ON s.arm = '%3$s' AND s.alpha = %4$s AND s.stratum = d.stratum
    WHERE d.u < CASE WHEN d.treatment = 1 THEN s.p_treated ELSE s.p_control END",
    stratum, seed, arm, format(alpha))
}
