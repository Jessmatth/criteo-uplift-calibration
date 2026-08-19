# ---- Build the confounded observational subsamples --------------------------
# The ad server in this study is a rule we write down, which is the only reason
# the residual bias at the end means anything. Construction:
#
#   1. Fit a baseline conversion model on the CONTROL arm only. Control users
#      never saw an ad, so their conversion probability is a clean measure of
#      pre-existing purchase intent, uncontaminated by the treatment effect.
#      This is the signal a real targeting system is built to chase.
#   2. Turn it into a targeting index g(X), deliberately nonlinear so that a
#      linear propensity model is genuinely misspecified rather than just noisy.
#   3. Cut users into K_STRATA equal-sized bins on g.
#   4. Retain treated users at a rate rising in the stratum and control users at
#      a rate falling in it. That is targeting-induced selection: the ad was
#      shown to the people who were going to buy anyway.
#
# Retention is constant WITHIN a stratum. That is the load-bearing choice: it
# leaves the within-stratum covariate distribution untouched by sampling, so the
# benchmark for each subsample is exactly the stratum-specific randomized
# effects reweighted to that subsample's stratum shares - no model involved.
source("R/00_config.R"); source("R/estimator.R"); source("R/confound.R")
library(DBI)

con <- db_connect()
set.seed(SEED)

# ---- 1. Baseline conversion model, control arm only -------------------------
message("Fitting baseline conversion model on the control arm ...")
# The sample must be taken INSIDE the control filter. DuckDB applies
# USING SAMPLE at the table scan, so sampling then filtering would draw from all
# 14M rows and leave only ~15% control - a mistake that fails silently, because
# the model still fits, just on a tenth of the intended data.
ctrl <- dbGetQuery(con, sprintf(
  "SELECT * FROM (SELECT %s, conversion FROM criteo WHERE treatment = 0)
   USING SAMPLE reservoir(1200000 ROWS) REPEATABLE (%d)",
  paste(FEATURES_ALL, collapse = ", "), SEED %% 1000))
fit <- glm(reformulate(FEATURES_ALL, "conversion"), data = ctrl,
           family = binomial())
b <- coef(fit)
message(sprintf("  n=%s, %d conversions, targeting model has %d coefficients",
                format(nrow(ctrl), big.mark = ","), sum(ctrl$conversion), length(b)))
rm(ctrl); invisible(gc())

# ---- 2/3. Targeting index and strata ---------------------------------------
# Strata are cut on the VALUE of g via cume_dist, not with ntile. This is not a
# stylistic choice. 6.68M distinct visible-feature combinations cover 13.98M
# rows, and the single most common combination accounts for 1,699,977 of them -
# 12% of the file share one identical visible profile. ntile splits such a tie
# group across several strata in scan order, and the scan order of this file is
# block-ordered by treatment (the first ~40% of rows are 100% treated). The
# result was a stratum that was 99.96% treated: retention would then have
# depended on file position rather than on X, which silently destroys both the
# selection-on-observables premise and the stratified benchmark. cume_dist gives
# every tied user the same stratum, so a stratum is a function of X alone.
# Strata are consequently unequal in size; §5 of protocol.md reports how uneven.
# Split the linear predictor into the part estimators can see and the part they
# cannot. In the "hidden" arm the second part drives targeting but is withheld
# from every model downstream, which is what makes the bias there irreducible.
lp_visible <- paste(sprintf("(%.10g) * %s", b[FEATURES_VISIBLE], FEATURES_VISIBLE),
                    collapse = " + ")
lp_hidden  <- sprintf("(%.10g) * %s", b[FEATURE_HIDDEN], FEATURE_HIDDEN)

message("Scoring all rows and cutting strata ...")
dbExecute(con, sprintf("
  CREATE OR REPLACE TABLE targeting AS
  WITH lp AS (
    SELECT user_id, %s AS lpv, %s AS lph, %s AS xa, %s AS xb FROM criteo
  ), z AS (
    SELECT user_id,
           (lpv - avg(lpv) OVER ()) / stddev_samp(lpv) OVER ()               AS z_observed,
           (lpv + %.10g * lph - avg(lpv + %.10g * lph) OVER ())
             / stddev_samp(lpv + %.10g * lph) OVER ()                        AS z_hidden,
           (xa - avg(xa) OVER ()) / stddev_samp(xa) OVER ()                  AS za,
           (xb - avg(xb) OVER ()) / stddev_samp(xb) OVER ()                  AS zb
    FROM lp
  ), g AS (
    -- Rounded to 9 decimals before binning. avg() and stddev_samp() over 14M
    -- rows are computed in parallel, and floating-point addition is not
    -- associative, so the standardization constants differ in the last bits
    -- between runs. The index is quadratic in z and therefore not monotone in
    -- it, so that jitter reorders values sitting exactly on a bin boundary and
    -- moves a few thousand rows between strata - enough to shift a headline
    -- number by several points. Rounding puts the run-to-run variation far
    -- below the bin width and makes the strata reproducible.
    SELECT user_id,
           round(z_observed + %.10g * (z_observed * z_observed - 1) + %.10g * za * zb, 9) AS g_observed,
           round(z_hidden   + %.10g * (z_hidden   * z_hidden   - 1) + %.10g * za * zb, 9) AS g_hidden
    FROM z
  )
  SELECT user_id, g_observed, g_hidden,
         least(%d, greatest(1, CAST(ceil(%d * cume_dist() OVER (ORDER BY g_observed)) AS INTEGER)))
           AS stratum_observed,
         least(%d, greatest(1, CAST(ceil(%d * cume_dist() OVER (ORDER BY g_hidden))   AS INTEGER)))
           AS stratum_hidden
  FROM g",
  lp_visible, lp_hidden, INTERACT_A, INTERACT_B,
  HIDDEN_LOADING, HIDDEN_LOADING, HIDDEN_LOADING,
  LAMBDA_SQ, KAPPA_INT, LAMBDA_SQ, KAPPA_INT,
  K_STRATA, K_STRATA, K_STRATA, K_STRATA))

# ---- 3b. Assert randomization survives stratification -----------------------
# A stratum is a function of X, and treatment was randomized, so the treated
# share must be flat across strata. This assertion is what catches a stratum
# that has picked up file ordering instead of covariate structure - the exact
# failure that ntile produced here, and one that is invisible downstream because
# every estimator still runs happily and returns a plausible-looking number.
for (arm in ARMS) {
  bal <- dbGetQuery(con, sprintf(
    "SELECT t.stratum_%1$s AS stratum, count(*) n, avg(c.treatment) trt_share,
            sum(CASE WHEN c.treatment=1 THEN 1 ELSE 0 END) n1,
            sum(CASE WHEN c.treatment=1 THEN c.conversion ELSE 0 END) x1_conversion,
            sum(CASE WHEN c.treatment=0 THEN 1 ELSE 0 END) n0,
            sum(CASE WHEN c.treatment=0 THEN c.conversion ELSE 0 END) x0_conversion
     FROM criteo c JOIN targeting t USING (user_id) GROUP BY 1 ORDER BY 1", arm))
  overall <- dbGetQuery(con, "SELECT avg(treatment) p FROM criteo")$p

  # Hard failure threshold. Structural contamination is gross - the ntile
  # version of this code produced strata ranging from 43% to 99.96% treated -
  # so 5 points is wide enough to tolerate the dataset's own mild imbalance and
  # narrow enough that no ordering artifact survives it.
  bad <- abs(bal$trt_share - overall) > 0.05
  if (any(bad)) {
    print(bal[bad, c("stratum", "n", "trt_share")])
    stop("treated share varies structurally across strata in arm '", arm,
         "': stratification is not a function of X alone")
  }
  # Milder deviations are reported, not tolerated silently. In this dataset the
  # top strata run ~2 points hot: users with high predicted purchase intent are
  # slightly likelier to sit in the treated arm. See the audit below.
  drift <- bal[abs(bal$trt_share - overall) > 0.01, ]
  message(sprintf("  arm %-8s %2d non-empty strata, largest %.1f%% of rows, treated share %.4f-%.4f%s",
                  arm, nrow(bal), 100 * max(bal$n) / sum(bal$n),
                  min(bal$trt_share), max(bal$trt_share),
                  if (nrow(drift)) sprintf("  [%d strata >1pt off]", nrow(drift)) else ""))
  write_table(bal, paste0("stratum_balance_", arm))
}

# ---- 3c. Does that imbalance move the benchmark? ----------------------------
# Criteo's holdout is not perfectly orthogonal to the features (max |SMD| across
# f0-f11 is 0.049, and the top targeting strata run ~2 points hot). If the
# randomization leaks, the crude difference in means is no longer guaranteed
# unbiased. The test: standardize the arm difference to the overall stratum
# distribution and see whether it moves. If crude and standardized agree, the
# leakage is too small to matter and the benchmark stands.
audit <- list()
for (arm in ARMS) {
  bal <- utils::read.csv(file.path(DIR_TABLES, paste0("stratum_balance_", arm, ".csv")))
  stopifnot(all(c("x1_conversion", "x0_conversion", "n1", "n0") %in% names(bal)))
  crude <- sum(bal$x1_conversion) / sum(bal$n1) - sum(bal$x0_conversion) / sum(bal$n0)
  std   <- stratified_ate(bal$n1, bal$x1_conversion, bal$n0, bal$x0_conversion, w = bal$n)
  audit[[arm]] <- data.frame(
    arm = arm,
    crude_itt = crude,
    standardized_itt = unname(std["estimate"]),
    standardized_se = unname(std["se"]),
    difference = unname(std["estimate"]) - crude,
    difference_in_ses = (unname(std["estimate"]) - crude) / unname(std["se"]))
}
audit <- do.call(rbind, audit)
# A crude ITT of exactly zero means a column was misread, not that the campaign
# had no effect: R/03 already measured it at 0.115pp on the same rows.
stopifnot(all(audit$crude_itt > 0), all(audit$standardized_itt > 0))
write_table(audit, "randomization_audit")
message("\n  Randomization audit (conversion, percentage points):")
for (i in seq_len(nrow(audit))) {
  message(sprintf("    %-8s crude ITT %.4f  vs  stratum-standardized %.4f  (%.2f SEs apart)",
                  audit$arm[i], 100 * audit$crude_itt[i],
                  100 * audit$standardized_itt[i], audit$difference_in_ses[i]))
}

# ---- 4. Stratum retention probabilities, solved per scenario ----------------
scen <- list(); truth <- list()
for (arm in ARMS) {
  sc <- dbGetQuery(con, sprintf("
    SELECT t.stratum_%1$s AS stratum,
           avg(t.g_%1$s) AS g_mean,
           sum(CASE WHEN c.treatment = 1 THEN 1 ELSE 0 END) AS n1,
           sum(CASE WHEN c.treatment = 1 THEN c.conversion ELSE 0 END) AS x1_conversion,
           sum(CASE WHEN c.treatment = 1 THEN c.visit ELSE 0 END)      AS x1_visit,
           sum(CASE WHEN c.treatment = 0 THEN 1 ELSE 0 END) AS n0,
           sum(CASE WHEN c.treatment = 0 THEN c.conversion ELSE 0 END) AS x0_conversion,
           sum(CASE WHEN c.treatment = 0 THEN c.visit ELSE 0 END)      AS x0_visit
    FROM criteo c JOIN targeting t USING (user_id)
    GROUP BY 1 ORDER BY 1", arm))
  m <- (sc$g_mean - mean(sc$g_mean)) / sd(sc$g_mean)

  for (alpha in ALPHA_GRID) {
    a1 <- solve_intercept(sc$n1, m, alpha, +1, N_OBS * TREATED_SHARE_TARGET)
    a0 <- solve_intercept(sc$n0, m, alpha, -1, N_OBS * (1 - TREATED_SHARE_TARGET))
    p1 <- retention(a1, m, alpha, +1)
    p0 <- retention(a0, m, alpha, -1)

    # The TRUE propensity score, known in closed form because we wrote the
    # selection rule. Within a stratum every unit shares one retention
    # probability per arm, so the probability of being treated GIVEN selection
    # is just the treated share of the retained units. Feeding this to IPTW and
    # AIPW gives an oracle reference: with the true score those estimators must
    # land on the benchmark. If they do not, the harness is wrong, not the data.
    e_true <- (sc$n1 * p1) / (sc$n1 * p1 + sc$n0 * p0)

    scen[[length(scen) + 1]] <- data.frame(
      arm = arm, alpha = alpha, stratum = sc$stratum,
      p_treated = p1, p_control = p0, e_true = e_true)

    # Ground truth for THIS subsample: stratum-specific randomized effects from
    # all 13.98M rows, reweighted to the subsample's stratum composition. ATE
    # weights use the subsample's total shares, ATT weights its treated shares.
    # ATE weights follow the subsample's whole stratum composition; ATT weights
    # follow its TREATED composition. Matching targets ATT, so grading it
    # against an ATE benchmark would charge it for a difference of estimand
    # rather than a failure of adjustment.
    w_ate <- (sc$n1 * p1 + sc$n0 * p0); w_ate <- w_ate / sum(w_ate)
    w_att <- (sc$n1 * p1);              w_att <- w_att / sum(w_att)

    for (o in c("conversion", "visit")) {
      x1 <- sc[[paste0("x1_", o)]]; x0 <- sc[[paste0("x0_", o)]]
      g_ate <- stratified_ate(sc$n1, x1, sc$n0, x0, w_ate)
      g_att <- stratified_ate(sc$n1, x1, sc$n0, x0, w_att)
      truth[[length(truth) + 1]] <- data.frame(
        scenario = scenario_id(arm, alpha), arm = arm, alpha = alpha,
        outcome = o, estimand = c("ATE", "ATT"),
        benchmark = c(g_ate["estimate"], g_att["estimate"]),
        benchmark_se = c(g_ate["se"], g_att["se"]),
        benchmark_lower = c(g_ate["lower"], g_att["lower"]),
        benchmark_upper = c(g_ate["upper"], g_att["upper"]),
        n_expected = sum(sc$n1 * p1 + sc$n0 * p0),
        treated_share = sum(sc$n1 * p1) / sum(sc$n1 * p1 + sc$n0 * p0),
        row.names = NULL)
    }
  }
}
scenarios <- do.call(rbind, scen)
dbWriteTable(con, "scenarios", scenarios, overwrite = TRUE)
benchmarks <- do.call(rbind, truth)
write_table(benchmarks, "scenario_benchmarks")
saveRDS(list(scenarios = scenarios, benchmarks = benchmarks, coef = b),
        "output/confound.rds")

# ---- Diagnostics: how strong is the induced selection? ----------------------
sel <- do.call(rbind, lapply(split(scenarios, list(scenarios$arm, scenarios$alpha)),
  function(d) data.frame(
    arm = d$arm[1], alpha = d$alpha[1],
    p_treated_lowest_stratum = d$p_treated[which.min(d$stratum)],
    p_treated_highest_stratum = d$p_treated[which.max(d$stratum)],
    treated_retention_ratio = d$p_treated[which.max(d$stratum)] /
                              d$p_treated[which.min(d$stratum)])))
sel <- sel[order(sel$arm, sel$alpha), ]
write_table(sel, "selection_strength")
print(sel, row.names = FALSE, digits = 3)

message("\n--- Benchmark ATE by scenario (conversion, percentage points) ---")
bt <- benchmarks[benchmarks$estimand == "ATE" & benchmarks$outcome == "conversion", ]
print(data.frame(scenario = bt$scenario,
                 benchmark_pp = round(100 * bt$benchmark, 5),
                 se_pp = round(100 * bt$benchmark_se, 5),
                 n = round(bt$n_expected)), row.names = FALSE)

# Explicit rather than on.exit: these scripts are source()d in sequence by
# run_all.R, where an on.exit handler fires on source()'s frame and leaves the
# next script holding an invalid connection.
dbDisconnect(con, shutdown = TRUE)
