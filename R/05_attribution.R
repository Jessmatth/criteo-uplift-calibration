# ---- Attribution credit versus incremental truth ----------------------------
# Everything here is computed from the real Criteo log. No journey is simulated
# and no touchpoint is invented; the synthetic multi-touch appendix lives in
# R/09_attribution_synthetic.R and is not used for any number in this file.
#
# The question a marketing team actually has to answer is not "what is the ATE".
# It is "how many of last month's conversions should be credited to this
# channel, and what did an incremental conversion cost". Last-touch answers it
# by crediting every conversion that followed an impression. The experiment
# answers it by differencing the arms. This file puts the two answers side by
# side, in conversions rather than percentage points.
source("R/00_config.R"); source("R/estimator.R")
library(DBI)

con <- db_connect(read_only = TRUE)

cells <- dbGetQuery(con, "
  SELECT treatment, exposure, count(*) n,
         sum(conversion) conversions, sum(visit) visits
  FROM criteo GROUP BY 1,2 ORDER BY 1,2")
trt_exp   <- cells[cells$treatment == 1 & cells$exposure == 1, ]
trt_unexp <- cells[cells$treatment == 1 & cells$exposure == 0, ]
ctl       <- cells[cells$treatment == 0, ]

# The audited, covariate-standardized ITT. R/03 showed the crude difference is
# inflated by residual imbalance in Criteo's holdout, so the standardized figure
# is the one a budget decision should use; the crude figure is carried alongside
# to show how much the audit is worth in dollars.
aud <- utils::read.csv(file.path(DIR_TABLES, "randomization_audit.csv"))
aud <- aud[aud$arm == "observed", ]

n_treated <- trt_exp$n + trt_unexp$n

rows <- list()
for (o in c("conversion", "visit")) {
  col <- if (o == "conversion") "conversions" else "visits"
  observed_total <- trt_exp[[col]] + trt_unexp[[col]]

  # Last touch: the ad was the last thing the user saw before converting, so the
  # ad gets the conversion. In a single-channel log that is every conversion
  # among exposed users.
  last_touch_credit <- trt_exp[[col]]

  # Truth: what the campaign added, from the randomized contrast.
  itt_crude <- observed_total / n_treated - ctl[[col]] / ctl$n
  itt_std   <- if (o == "conversion") aud$standardized_itt else NA_real_
  incr_crude <- itt_crude * n_treated
  incr_std   <- itt_std * n_treated

  rows[[o]] <- data.frame(
    outcome = o,
    observed_events_treated = observed_total,
    last_touch_credit = last_touch_credit,
    last_touch_share_of_events = last_touch_credit / observed_total,
    incremental_crude = incr_crude,
    incremental_audited = incr_std,
    true_share_of_events_crude = incr_crude / observed_total,
    true_share_of_events_audited = incr_std / observed_total,
    overstatement_vs_crude = last_touch_credit / incr_crude,
    overstatement_vs_audited = last_touch_credit / incr_std)
}
credit <- do.call(rbind, rows)
write_table(credit, "attribution_credit")

# ---- Estimand confusion, priced ---------------------------------------------
# Three numbers get called "the effect of the ad". They differ by 45x. Which one
# a team reports decides whether the campaign looks profitable.
b <- readRDS("output/benchmark.rds")$benchmark
conv <- b[b$outcome == "conversion", ]
cmp <- data.frame(
  reported_as = c("last-touch / exposed vs unexposed",
                  "effect of seeing an ad (CACE, IV)",
                  "effect of running the campaign (ITT)"),
  per_user_effect_pp = 100 * c(conv$estimate[conv$estimand == "NAIVE exposed vs unexposed"],
                               conv$estimate[conv$estimand == "CACE (IV on exposure)"],
                               conv$estimate[conv$estimand == "ITT (randomized)"]),
  applies_to = c("the 3.6% of treated users who saw an ad",
                 "the users whose exposure the campaign changed",
                 "every user the campaign could reach"),
  is_causal = c("no", "yes", "yes"))
write_table(cmp, "estimand_comparison")

cc <- credit["conversion", ]
message("\n--- Touchpoint credit, conversions (real data only) ---")
message(sprintf("  Conversions among ad-eligible users        %s", format(round(cc$observed_events_treated), big.mark = ",")))
message(sprintf("  Credited to the ad by last-touch           %s  (%.1f%% of them)",
                format(round(cc$last_touch_credit), big.mark = ","), 100 * cc$last_touch_share_of_events))
message(sprintf("  Actually incremental (crude ITT)           %s  (%.1f%%)",
                format(round(cc$incremental_crude), big.mark = ","), 100 * cc$true_share_of_events_crude))
message(sprintf("  Actually incremental (audited ITT)         %s  (%.1f%%)",
                format(round(cc$incremental_audited), big.mark = ","), 100 * cc$true_share_of_events_audited))
message(sprintf("\n  Last-touch overstates the ad's contribution by %.2fx (crude) / %.2fx (audited).",
                cc$overstatement_vs_crude, cc$overstatement_vs_audited))
message(sprintf("  Cost per incremental conversion is therefore %.2fx what last-touch reports.",
                cc$overstatement_vs_audited))

# Explicit rather than on.exit: these scripts are source()d in sequence by
# run_all.R, where an on.exit handler fires on source()'s frame and leaves the
# next script holding an invalid connection.
dbDisconnect(con, shutdown = TRUE)
