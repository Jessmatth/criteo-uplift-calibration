#!/usr/bin/env Rscript
# Full pipeline from the repository root:  Rscript run_all.R
#
# Tunables:
#   N_OBS=500000 Rscript run_all.R    # smaller confounded subsamples, faster
#   N_FOLDS=2    Rscript run_all.R    # fewer cross-fitting folds
#
# First run downloads a 297 MB archive into data/raw/ and builds a DuckDB
# database from it. At the default N_OBS=2000000 the whole pipeline takes about
# 26 minutes on 8 cores after the download, of which R/06_estimators.R is 25.

steps <- c(
  "R/01_download.R",             # fetch and verify the Criteo archive
  "R/02_build_analytic.R",       # stage into DuckDB, validate the experiment
  "R/03_benchmark.R",            # randomized ground truth, and audit it
  "R/04_confound.R",             # build the confounded observational subsamples
  "R/05_attribution.R",          # last-touch credit vs incremental truth
  "R/06_estimators.R",           # IPTW, matching, AIPW across the sweep
  "R/07_calibrate.R",            # grade every estimate against the benchmark
  "R/08_guardrails.R",           # SLA pass/fail scorecard
  "R/09_attribution_synthetic.R", # APPENDIX: simulated multi-touch journeys
  "R/10_figures.R"
)

for (.step in steps) {
  message("\n", strrep("=", 72), "\n  ", .step, "\n", strrep("=", 72))
  # Dotted name on purpose: source() evaluates each script in the global
  # environment, so a plain `t0` here is fair game for any script to overwrite.
  .step_start <- Sys.time()
  source(.step, echo = FALSE)
  message(sprintf("  [%.1fs]",
                  as.numeric(difftime(Sys.time(), .step_start, units = "secs"))))
}

message("\nRunning tests ...")
testthat::test_dir("tests/testthat", stop_on_failure = TRUE)
message("\nPipeline complete. Tables in output/tables/, figures in output/figures/.")
