# ---- Experimental ground truth ----------------------------------------------
# Three numbers that the rest of the repository is measured against, all
# computed on the full randomized sample:
#   ITT   - effect of being made ad-eligible, on everyone. What a holdout test
#           actually measures, and what a campaign budget decision turns on.
#   CACE  - effect of actually seeing an ad, among the users whose exposure the
#           campaign changed. The estimand attribution systems believe they are
#           reporting.
#   NAIVE - exposed users versus unexposed users. Not a causal quantity at all.
#           It is here because it is the comparison the industry ships.
source("R/00_config.R"); source("R/estimator.R")
library(DBI)

con <- db_connect(read_only = TRUE)

b <- dbGetQuery(con, read_sql("sql/02_benchmark.sql"))
arm1 <- b[b$treatment == 1, ]; arm0 <- b[b$treatment == 0, ]

rows <- list()
for (o in c("conversion", "visit")) {
  itt <- newcombe_rd(arm1[[paste0(o, "s")]], arm1$n, arm0[[paste0(o, "s")]], arm0$n)
  cace <- cace_wald_moments(
    n1 = arm1$n, ybar1 = arm1[[paste0("ybar_", o)]], dbar1 = arm1$dbar_exposure,
    ydbar1 = arm1[[paste0("ydbar_", o)]],
    n0 = arm0$n, ybar0 = arm0[[paste0("ybar_", o)]], dbar0 = arm0$dbar_exposure,
    ydbar0 = arm0[[paste0("ydbar_", o)]])

  rows[[length(rows) + 1]] <- data.frame(
    outcome = o, estimand = "ITT (randomized)",
    estimate = unname(itt["estimate"]),
    lower = unname(itt["lower.lower"]), upper = unname(itt["upper.upper"]),
    interpretation = "effect of ad eligibility on all users")
  rows[[length(rows) + 1]] <- data.frame(
    outcome = o, estimand = "CACE (IV on exposure)",
    estimate = unname(cace["estimate"]),
    lower = unname(cace["lower"]), upper = unname(cace["upper"]),
    interpretation = "effect of actually seeing an ad, among compliers")
}

# The comparison an attribution system makes: converted-and-exposed against
# everyone else in the treated population. No randomization protects it.
cells <- dbGetQuery(con, "
  SELECT exposure, count(*) n, sum(conversion) conversions, sum(visit) visits
  FROM criteo WHERE treatment = 1 GROUP BY exposure ORDER BY exposure")
ex <- cells[cells$exposure == 1, ]; un <- cells[cells$exposure == 0, ]
for (o in c("conversion", "visit")) {
  col <- paste0(o, "s")
  nv <- newcombe_rd(ex[[col]], ex$n, un[[col]], un$n)
  rows[[length(rows) + 1]] <- data.frame(
    outcome = o, estimand = "NAIVE exposed vs unexposed",
    estimate = unname(nv["estimate"]),
    lower = unname(nv["lower.lower"]), upper = unname(nv["upper.upper"]),
    interpretation = "not causal: the ad server chose who was exposed")
}
bench <- do.call(rbind, rows)
write_table(bench, "benchmark_randomized")

# ---- What the selection looks like, measured not assumed ---------------------
compliance <- arm1$dbar_exposure
cb <- complier_baseline(p_control = arm0$ybar_conversion,
                        p_treated_unexposed = un$conversions / un$n,
                        compliance = compliance)
sel <- data.frame(
  quantity = c("compliance rate (share of treated actually exposed)",
               "conversion rate, treated & unexposed (never-takers)",
               "conversion rate, control arm (mixture)",
               "implied untreated conversion rate of compliers",
               "complier baseline / never-taker baseline"),
  value = c(compliance, un$conversions / un$n, arm0$ybar_conversion, cb,
            cb / (un$conversions / un$n)))
write_table(sel, "selection_evidence")

# ---- Audit the experiment before trusting it --------------------------------
# A randomized benchmark is only a benchmark if the randomization held. Under
# clean assignment, standardizing the arm difference to the marginal
# distribution of ANY covariate must leave it unchanged, because the covariate
# is balanced by construction. This check needs no model: bin each raw feature,
# recompute the arm difference within bins, and reweight to the pooled bin
# sizes. Any systematic movement is evidence of imbalance.
aud <- do.call(rbind, lapply(FEATURES_ALL, function(f) {
  d <- dbGetQuery(con, sprintf("
    SELECT b, count(*) n,
           sum(CASE WHEN treatment=1 THEN 1 ELSE 0 END) n1,
           sum(CASE WHEN treatment=1 THEN conversion ELSE 0 END) x1,
           sum(CASE WHEN treatment=0 THEN 1 ELSE 0 END) n0,
           sum(CASE WHEN treatment=0 THEN conversion ELSE 0 END) x0
    FROM (SELECT treatment, conversion,
                 least(20, greatest(1, CAST(ceil(20 * cume_dist() OVER (ORDER BY %s))
                       AS INTEGER))) AS b
          FROM criteo) GROUP BY 1 ORDER BY 1", f))
  d <- d[d$n1 > 0 & d$n0 > 0, ]
  crude <- sum(d$x1) / sum(d$n1) - sum(d$x0) / sum(d$n0)
  st <- stratified_ate(d$n1, d$x1, d$n0, d$x0, w = d$n)
  data.frame(feature = f, bins = nrow(d), crude_itt = crude,
             standardized_itt = unname(st["estimate"]),
             ses_apart = (unname(st["estimate"]) - crude) / unname(st["se"]))
}))
write_table(aud, "randomization_audit_marginal")

worst <- aud[which.min(aud$ses_apart), ]
message(sprintf("\n--- Randomization audit ---"))
message(sprintf("  Standardizing on a single raw feature moves the ITT from %.4f%% to as low as %.4f%%",
                100 * worst$crude_itt, 100 * min(aud$standardized_itt)))
message(sprintf("  Largest shift: %s, %.1f SEs. %d of %d features move it DOWN.",
                worst$feature, worst$ses_apart, sum(aud$ses_apart < 0), nrow(aud)))
message("  Criteo's holdout is not orthogonal to the features. The crude ITT is")
message("  inflated; R/04_confound.R standardizes it, and every scenario benchmark")
message("  in this study is covariate-standardized rather than crude.")

saveRDS(list(benchmark = bench, selection = sel, arms = b),
        "output/benchmark.rds")

fmt <- function(x) sprintf("%.4f%%", 100 * x)
message("\n--- Randomized ground truth (conversion) ---")
for (i in which(bench$outcome == "conversion")) {
  message(sprintf("  %-28s %s  [%s, %s]", bench$estimand[i],
                  fmt(bench$estimate[i]), fmt(bench$lower[i]), fmt(bench$upper[i])))
}
message(sprintf("\n  Compliers' untreated conversion rate is %.1fx the never-takers'.",
                cb / (un$conversions / un$n)))
message("  That ratio is the targeting selection, read directly off the experiment.")

# Explicit rather than on.exit: these scripts are source()d in sequence by
# run_all.R, where an on.exit handler fires on source()'s frame and leaves the
# next script holding an invalid connection.
dbDisconnect(con, shutdown = TRUE)
