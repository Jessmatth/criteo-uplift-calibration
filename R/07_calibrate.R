# ---- Calibrate every estimate against the experiment ------------------------
# The single question this file answers, for each estimator and each severity of
# targeting: how much of the gap between a naive readout and the truth does it
# actually close, and does its stated uncertainty tell the truth about that.
source("R/00_config.R"); source("R/estimator.R")

est   <- utils::read.csv(file.path(DIR_TABLES, "estimates_raw.csv"))
bench <- utils::read.csv(file.path(DIR_TABLES, "scenario_benchmarks.csv"))
diag  <- utils::read.csv(file.path(DIR_TABLES, "estimator_diagnostics.csv"))

cal <- merge(est, bench[, c("scenario", "outcome", "estimand", "benchmark",
                            "benchmark_se", "benchmark_lower", "benchmark_upper")],
             by = c("scenario", "outcome", "estimand"), all.x = TRUE)
stopifnot(!any(is.na(cal$benchmark)))

cal$bias          <- cal$estimate - cal$benchmark
cal$relative_bias <- cal$bias / cal$benchmark
# Is the observed bias larger than the noise in measuring it? Estimator and
# benchmark are computed on overlapping data, so treating them as independent
# overstates the variance of the difference - this is deliberately conservative.
cal$bias_se       <- sqrt(cal$se^2 + cal$benchmark_se^2)
cal$bias_z        <- cal$bias / cal$bias_se
cal$ci_covers_benchmark <- cal$lower <= cal$benchmark & cal$benchmark <= cal$upper
cal$ci_width_ratio <- (cal$upper - cal$lower) /
                      (cal$benchmark_upper - cal$benchmark_lower)

# Share of the naive gap removed. 1 means the estimator recovered the
# experiment; 0 means it did no better than the unadjusted difference; negative
# means adjustment moved the answer further from the truth, which does happen.
naive <- cal[cal$estimator == "naive (unadjusted)",
             c("scenario", "outcome", "estimand", "bias")]
names(naive)[4] <- "naive_bias"
cal <- merge(cal, naive, by = c("scenario", "outcome", "estimand"), all.x = TRUE)
cal$gap_closed <- 1 - abs(cal$bias) / abs(cal$naive_bias)

cal <- merge(cal, diag[, c("scenario", "ps_spec", "max_abs_smd_weighted",
                           "extreme_ps_share", "ess_ratio")],
             by = c("scenario", "ps_spec"), all.x = TRUE)
cal <- cal[order(cal$outcome, cal$arm, cal$alpha, cal$estimator, cal$ps_spec,
                 cal$om_spec), ]
write_table(cal, "calibration")
saveRDS(cal, "output/calibration.rds")

# ---- Headline view: conversion, ATE estimators ------------------------------
# Relative bias leads, not share-of-gap-closed. The naive readout is 350-590%
# off, so closing 90% of that gap still leaves a 50% overstatement - a metric
# that reads as a triumph while the answer is still unusable.
show <- cal[cal$outcome == "conversion" & cal$estimand == "ATE" & cal$alpha > 0, ]
show$label <- ifelse(is.na(show$ps_spec) & is.na(show$om_spec), show$estimator,
              ifelse(is.na(show$om_spec), paste0(show$estimator, " [ps=", show$ps_spec, "]"),
              ifelse(is.na(show$ps_spec), paste0(show$estimator, " [om=", show$om_spec, "]"),
                     paste0(show$estimator, " [ps=", show$ps_spec, ", om=", show$om_spec, "]"))))

wide <- function(d, col) {
  w <- reshape(d[, c("label", "alpha", col)], idvar = "label", timevar = "alpha",
               direction = "wide")
  names(w) <- sub(paste0(col, "\\."), "alpha=", names(w))
  w[-1] <- round(100 * w[-1], 0)
  w[order(w$label), ]
}
for (a in ARMS) {
  message(sprintf("\n=== %s arm: relative bias %%, conversion (0 = recovered the experiment) ===",
                  toupper(a)))
  print(wide(show[show$arm == a, ], "relative_bias"), row.names = FALSE)
}

# ---- What a practitioner could actually achieve ------------------------------
# The oracle propensity score is known only because this study wrote the
# selection rule. Excluding it, the best any specification manages is the
# honest ceiling on observational ad measurement here.
feas <- show[is.na(show$ps_spec) | show$ps_spec != "oracle", ]
feas <- feas[feas$estimator != "naive (unadjusted)", ]
best <- do.call(rbind, lapply(split(feas, list(feas$arm, feas$alpha), drop = TRUE),
  function(d) {
    i <- which.min(abs(d$relative_bias))
    nv <- show$relative_bias[show$arm == d$arm[1] & show$alpha == d$alpha[1] &
                             show$estimator == "naive (unadjusted)"]
    data.frame(arm = d$arm[1], alpha = d$alpha[1],
               naive_relative_bias = nv,
               best_feasible_spec = d$label[i],
               best_relative_bias = d$relative_bias[i],
               ci_covers_benchmark = d$ci_covers_benchmark[i])
  }))
best <- best[order(best$arm, best$alpha), ]
write_table(best, "best_feasible")
message("\n=== Best feasible specification (oracle excluded) ===")
print(data.frame(arm = best$arm, alpha = best$alpha,
                 naive_bias_pct = round(100 * best$naive_relative_bias, 0),
                 best_bias_pct = round(100 * best$best_relative_bias, 0),
                 spec = best$best_feasible_spec,
                 covers = best$ci_covers_benchmark), row.names = FALSE)
