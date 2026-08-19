# ---- Pass/fail harness ------------------------------------------------------
# Turns the calibration table into a decision. Thresholds come from
# GUARDRAILS in R/00_config.R and were fixed before any estimate existed;
# changing one after seeing results is an amendment, and protocol.md says so.
#
# A FAIL here is usually a finding, not a defect. The interesting output is not
# which estimators pass but WHERE they stop passing: how much targeting
# pressure, and how much unobserved confounding, an observational readout can
# absorb before it should no longer inform a budget.
source("R/00_config.R"); source("R/estimator.R"); source("R/guardrails.R")

cal <- readRDS("output/calibration.rds")
diag <- utils::read.csv(file.path(DIR_TABLES, "estimator_diagnostics.csv"))

score <- do.call(rbind, lapply(seq_len(nrow(cal)), function(i) {
  r <- cal[i, ]
  d <- if (!is.na(r$ps_spec))
         diag[diag$scenario == r$scenario & diag$ps_spec == r$ps_spec, ] else NULL
  cbind(r[, c("scenario", "arm", "alpha", "outcome", "estimator", "estimand",
              "ps_spec", "om_spec", "estimate", "lower", "upper", "benchmark",
              "gap_closed")],
        evaluate_row(r, r, d, GUARDRAILS))
}))
write_table(score, "guardrail_scorecard")
saveRDS(score, "output/scorecard.rds")

# ---- Harness self-check -----------------------------------------------------
# With the true propensity score and no induced confounding, IPTW must land on
# the benchmark. This is a sharp test rather than a hopeful one: Hajek IPTW
# weighted by the stratum-level true score is algebraically the stratum-
# standardized estimator, which is exactly how the benchmark is computed. If
# these two disagree, the construction is wrong and nothing else can be read.
#
# The test deliberately uses IPTW and not AIPW. AIPW's augmentation term adjusts
# at FEATURE resolution via an outcome model of the real conversion surface,
# whose functional form is unknown, while the benchmark and the oracle score are
# both defined at STRATUM resolution. Where the underlying randomization is
# imperfect - and here it is, see the audit in R/03 - those two resolutions
# answer slightly different questions, so an AIPW discrepancy at alpha = 0
# indicts the outcome model, not the construction. That discrepancy is itself
# reported below rather than hidden.
selfcheck <- score[score$alpha == 0 & score$ps_spec == "oracle" &
                   score$estimator == "IPTW (stabilized)" &
                   score$estimand == "ATE" & score$outcome == "conversion", ]
message("\n--- Harness self-check (alpha=0, true propensity score, IPTW) ---")
for (i in seq_len(nrow(selfcheck))) {
  message(sprintf("  %-14s estimate %.5f%%  benchmark %.5f%%  rel.bias %+.1f%%  covers: %s",
                  selfcheck$scenario[i], 100 * selfcheck$estimate[i],
                  100 * selfcheck$benchmark[i], 100 * selfcheck$relative_bias[i],
                  selfcheck$ci_covers[i]))
}
if (any(!selfcheck$ci_covers)) {
  stop("harness self-check failed: with the true propensity score and no induced ",
       "confounding, IPTW does not cover the benchmark. The construction is wrong, ",
       "not the estimator.")
}

# What the naive readout and AIPW do in that same no-confounding cell. The only
# thing there is to correct at alpha = 0 is Criteo's own residual imbalance, so
# this is a clean measure of how much of it each approach removes.
keep <- score$alpha == 0 & score$estimand == "ATE" &
        score$outcome == "conversion" &
        (score$estimator == "naive (unadjusted)" |
         (!is.na(score$ps_spec) & score$ps_spec == "oracle" &
          !is.na(score$om_spec) & score$om_spec == "rich"))
zero <- score[which(keep), ]
message("\n  With no induced confounding, correcting only Criteo's own imbalance:")
for (i in seq_len(nrow(zero))) {
  message(sprintf("    %-14s %-28s rel.bias %+6.1f%%  covers: %s",
                  zero$scenario[i], zero$estimator[i],
                  100 * zero$relative_bias[i], zero$ci_covers[i]))
}

# ---- Summary ----------------------------------------------------------------
sc <- score[score$outcome == "conversion", ]
sc$pass <- sc$verdict == "PASS"
summ <- aggregate(pass ~ estimator + arm + alpha, data = sc, FUN = mean)
summ$pass <- round(100 * summ$pass, 0)
w <- reshape(summ, idvar = c("estimator", "arm"), timevar = "alpha",
             direction = "wide")
names(w) <- sub("pass\\.", "alpha=", names(w))
message("\n--- % of specifications passing all guardrails (conversion) ---")
print(w[order(w$arm, w$estimator), ], row.names = FALSE)

# ---- Does a commoner outcome change the verdict? ----------------------------
# Conversion runs at 0.3%; visit at 4.85%. If specifications fail on conversion
# but pass on visit, the binding constraint was precision, not bias - a
# different problem with a different fix (more data, not better adjustment).
byout <- do.call(rbind, lapply(split(score, list(score$outcome, score$alpha), drop = TRUE),
  function(d) data.frame(outcome = d$outcome[1], alpha = d$alpha[1],
                         pass_rate = mean(d$verdict == "PASS"),
                         width_failures = mean(!d$ci_width_pass, na.rm = TRUE),
                         bias_failures = mean(!d$relative_bias_pass, na.rm = TRUE))))
byout <- byout[order(byout$outcome, byout$alpha), ]
write_table(byout, "guardrails_by_outcome")
message("\n--- Pass rate and failure mode by outcome ---")
print(data.frame(outcome = byout$outcome, alpha = byout$alpha,
                 pass_pct = round(100 * byout$pass_rate),
                 ci_width_fail_pct = round(100 * byout$width_failures),
                 bias_fail_pct = round(100 * byout$bias_failures)), row.names = FALSE)

fail_reasons <- score[score$verdict == "FAIL" & score$outcome == "conversion", ]
tally <- data.frame(
  check = c("relative bias", "CI covers benchmark", "CI width", "positivity",
            "balance", "ESS"),
  failures = c(sum(!fail_reasons$relative_bias_pass, na.rm = TRUE),
               sum(!fail_reasons$ci_covers_pass, na.rm = TRUE),
               sum(!fail_reasons$ci_width_pass, na.rm = TRUE),
               sum(!fail_reasons$positivity_pass, na.rm = TRUE),
               sum(!fail_reasons$balance_pass, na.rm = TRUE),
               sum(!fail_reasons$ess_pass, na.rm = TRUE)))
write_table(tally, "guardrail_failure_reasons")
message("\n--- Why specifications fail (conversion) ---")
print(tally, row.names = FALSE)
