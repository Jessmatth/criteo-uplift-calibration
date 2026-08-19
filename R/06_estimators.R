# ---- Run every estimator on every confounded subsample ----------------------
# Four families, on identical data, graded later against the same benchmark:
#
#   naive   - difference in means. What an unadjusted dashboard reports.
#   IPTW    - reweight by the inverse propensity of exposure. Right only if the
#             propensity model is right.
#   matching- 1:1 on the linear propensity score, with replacement. Targets ATT,
#             not ATE, and is graded against the ATT benchmark accordingly.
#   AIPW    - doubly robust. Right if EITHER the propensity model or the outcome
#             model is right.
#
# Each is run with a deliberately misspecified propensity model (linear in the
# features, when selection is quadratic plus an interaction), with the correct
# functional form, and with the TRUE propensity score known from construction.
# Crossing those against a parsimonious and a rich outcome model is what turns
# "AIPW is doubly robust" from a claim into a measurement.
source("R/00_config.R"); source("R/estimator.R"); source("R/confound.R")
library(DBI)

con <- db_connect(read_only = TRUE)
coefs <- readRDS("output/confound.rds")$coef
OUTCOMES <- c("conversion", "visit")

results <- list(); diags <- list()

for (arm in ARMS) for (alpha in ALPHA_GRID) {
  sid <- scenario_id(arm, alpha)
  t0 <- Sys.time()
  d <- dbGetQuery(con, obs_sample_sql(arm, alpha, SEED, K_STRATA))
  n <- nrow(d); treat <- d$treatment

  # ---- Design matrices ------------------------------------------------------
  # Estimators see FEATURES_VISIBLE only, in BOTH arms. In the hidden arm the
  # targeting index also loads on f0, which is withheld here; that omission is
  # the entire difference between the two arms.
  V <- d[, FEATURES_VISIBLE, drop = FALSE]
  z  <- as.numeric(scale(as.matrix(V) %*% coefs[FEATURES_VISIBLE]))
  za <- as.numeric(scale(d[[INTERACT_A]])); zb <- as.numeric(scale(d[[INTERACT_B]]))
  X <- V; X$.z2 <- z^2; X$.zab <- za * zb

  f_lin  <- reformulate(FEATURES_VISIBLE, ".treat")
  f_flex <- reformulate(c(FEATURES_VISIBLE, ".z2", ".zab"), ".treat")
  o_lin  <- reformulate(FEATURES_VISIBLE, ".y")
  o_rich <- reformulate(c(FEATURES_VISIBLE, ".z2", ".zab"), ".y")

  folds <- make_folds(n, N_FOLDS, SEED + as.integer(alpha * 7) + nchar(arm))

  # ---- Propensity scores ----------------------------------------------------
  # `oracle` is not available to a real analyst. It is here as a control: with
  # the true score, IPTW and AIPW must recover the benchmark, so a failure there
  # indicts this harness rather than the estimator.
  ps_raw <- list(
    linear   = cross_fit_ps(X, treat, f_lin,  folds),
    flexible = cross_fit_ps(X, treat, f_flex, folds),
    oracle   = d$e_true)
  # Truncation is applied for estimation; the untruncated scores are kept so the
  # positivity diagnostic reports how much mass truncation actually moved,
  # rather than reporting on the already-clipped values and finding nothing.
  ps <- lapply(ps_raw, truncate_ps, bounds = PS_TRUNCATE)

  # ---- Diagnostics ----------------------------------------------------------
  smd_unadj <- max(abs(sapply(FEATURES_VISIBLE, function(f) smd(d[[f]], treat))))
  for (nm in names(ps)) {
    w <- stabilized_weights(treat, ps[[nm]])
    diags[[length(diags) + 1]] <- data.frame(
      scenario = sid, arm = arm, alpha = alpha, ps_spec = nm,
      max_abs_smd_unadjusted = smd_unadj,
      max_abs_smd_weighted = max(abs(sapply(FEATURES_VISIBLE,
                                            function(f) smd(d[[f]], treat, w)))),
      extreme_ps_share = mean(ps_raw[[nm]] < PS_TRUNCATE[1] | ps_raw[[nm]] > PS_TRUNCATE[2]),
      ess_ratio = ess(w) / n,
      min_ps = min(ps_raw[[nm]]), max_ps = max(ps_raw[[nm]]))
  }

  add <- function(outcome, estimator, estimand, r, ps_spec = NA, om_spec = NA) {
    results[[length(results) + 1]] <<- data.frame(
      scenario = sid, arm = arm, alpha = alpha, outcome = outcome,
      estimator = estimator, estimand = estimand,
      ps_spec = ps_spec, om_spec = om_spec,
      estimate = unname(r["estimate"]), se = unname(r["se"]),
      lower = unname(r["lower"]), upper = unname(r["upper"]))
  }

  for (o in OUTCOMES) {
    y <- d[[o]]

    # naive
    nv <- newcombe_rd(sum(y[treat == 1]), sum(treat == 1),
                      sum(y[treat == 0]), sum(treat == 0))
    add(o, "naive (unadjusted)", "ATE",
        c(estimate = unname(nv["estimate"]),
          se = (unname(nv["upper.upper"]) - unname(nv["lower.lower"])) / (2 * 1.96),
          lower = unname(nv["lower.lower"]), upper = unname(nv["upper.upper"])))

    # IPTW and matching, per propensity spec
    for (nm in names(ps)) {
      add(o, "IPTW (stabilized)", "ATE", iptw_ate(y, treat, ps[[nm]]), ps_spec = nm)
      lps <- stats::qlogis(ps[[nm]])
      m <- match_1to1(lps, treat, MATCH_CALIPER_SD)
      add(o, "PS matching 1:1", "ATT", matched_att(y, m), ps_spec = nm)
    }

    # outcome models, then g-computation and the full AIPW grid
    om <- list(parsimonious = cross_fit_om(X, y, treat, o_lin,  folds),
               rich         = cross_fit_om(X, y, treat, o_rich, folds))
    for (onm in names(om)) {
      add(o, "g-computation", "ATE", gcomp_ate(om[[onm]]$mu1, om[[onm]]$mu0),
          om_spec = onm)
      for (nm in names(ps)) {
        add(o, "AIPW (cross-fit)", "ATE",
            aipw_ate(y, treat, ps[[nm]], om[[onm]]$mu1, om[[onm]]$mu0),
            ps_spec = nm, om_spec = onm)
      }
    }
    rm(om)
  }
  rm(d, V, X, ps, ps_raw); invisible(gc())
  message(sprintf("  %-14s n=%s  [%.0fs]", sid, format(n, big.mark = ","),
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

estimates <- do.call(rbind, results)
diagnostics <- do.call(rbind, diags)
write_table(estimates, "estimates_raw")
write_table(diagnostics, "estimator_diagnostics")
saveRDS(list(estimates = estimates, diagnostics = diagnostics), "output/estimates.rds")
message(sprintf("\n%d estimates across %d scenarios.", nrow(estimates),
                length(unique(estimates$scenario))))

# Explicit rather than on.exit: these scripts are source()d in sequence by
# run_all.R, where an on.exit handler fires on source()'s frame and leaves the
# next script holding an invalid connection.
dbDisconnect(con, shutdown = TRUE)
