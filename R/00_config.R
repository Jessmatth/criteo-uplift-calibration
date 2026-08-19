# ---- Project configuration --------------------------------------------------
# Criteo Uplift Prediction dataset v2.1 (Criteo AI Lab, AdKDD 2018).
#
# The file is the log of a real randomized incrementality test: a random ~15%
# of users were withheld from Criteo's ad targeting entirely, the remaining
# ~85% were eligible to be served ads. Because assignment was randomized, the
# incremental effect of the campaign is identified without assumptions. That
# number is the benchmark every observational estimate in this repository is
# graded against.
#
# NOTE ON SIZE: v2.1 contains 13,979,592 rows. The frequently cited "25M rows"
# refers to the original 2018 v2.0 release, which Criteo replaced with this
# subsampled version. Verified locally by row count, not taken from the paper.

DIR_RAW     <- "data/raw"
DIR_DERIVED <- "data/derived"
DIR_TABLES  <- "output/tables"
DIR_FIGURES <- "output/figures"
PATH_DUCKDB <- file.path(DIR_DERIVED, "criteo.duckdb")
PATH_GZ     <- file.path(DIR_RAW, "criteo-uplift-v2.1.csv.gz")

# Criteo serves this through a redirector; the 302 lands on Azure blob storage.
DATA_URL <- "http://go.criteo.net/criteo-research-uplift-v2.1.csv.gz"
DATA_MD5 <- "d2236769ef69e9be52556110102911ec"
DATA_BYTES <- 311422618
DATA_ROWS  <- 13979592

SEED <- 20260819

# ---- Column roles -----------------------------------------------------------
# f0-f11 are anonymized dense continuous features; Criteo never published their
# semantics. f0 is reserved as the HIDDEN confounder: in the "hidden" arm it
# drives ad targeting but is withheld from every estimator, which is what makes
# the residual bias in that arm irreducible rather than a specification failure.
FEATURE_HIDDEN   <- "f0"
FEATURES_VISIBLE <- paste0("f", 1:11)
FEATURES_ALL     <- c(FEATURE_HIDDEN, FEATURES_VISIBLE)

# The analyzed treatment is the RANDOMIZED assignment flag, not `exposure`.
# `exposure` is post-randomization (only treated users can be exposed), so it is
# never used as a treatment variable in the calibration study. It gets its own
# analysis in R/05_attribution.R, where it is the thing being critiqued.
COL_TREATMENT <- "treatment"
COL_OUTCOME   <- "conversion"
COL_OUTCOME_2 <- "visit"

# ---- Confounding construction ----------------------------------------------
# Users are ranked by a targeting index and cut into K_STRATA equal-sized bins.
# Retention into the observational sample is a function of the STRATUM, not of
# the continuous index. That choice is load-bearing: because every unit in a
# stratum shares one retention probability, the within-stratum covariate
# distribution survives sampling untouched, so the stratified randomized
# estimate of the ground truth is exactly unbiased rather than approximately so.
# See protocol.md section 5.
K_STRATA <- 50L

# Selection strength. The retention ratio between the most and least targeted
# stratum runs roughly exp(2.9 * alpha), so this grid spans no selection, mild
# (~2x), moderate (~4x) and severe (~18x) targeting.
#
# The grid was rescaled after a first pass: at alpha = 1 the true propensity
# score reached 0.9996 and 11.6% of units fell outside the truncation bounds,
# i.e. whole strata had almost no controls left. That is a positivity violation,
# and under one no estimator can recover the effect - the guardrails correctly
# failed every cell, but a sweep in which every non-zero cell is infeasible
# measures nothing. alpha = 1 is retained as the severe end precisely so the
# breakdown is visible rather than assumed.
#
# alpha = 0 is a negative control: no confounding is induced, so every estimator
# including the naive difference must land on the benchmark. If it does not, the
# harness is broken, not the estimator.
ALPHA_GRID <- c(0, 0.25, 0.5, 1)

# Two targeting regimes.
#   observed - the index is built from FEATURES_VISIBLE only. A correctly
#              specified estimator has everything it needs; bias is recoverable.
#   hidden   - the index also loads on FEATURE_HIDDEN, which no estimator sees.
#              No amount of adjustment closes this gap. This is the arm that
#              answers "when should I not trust an observational readout?"
ARMS <- c("observed", "hidden")
# Weight on the hidden feature's contribution in the "hidden" arm, chosen by
# how much of the targeting index it makes unreconstructable:
#
#   loading    corr(g_obs, g_hidden)    share of index not recoverable from f1-f11
#         1                    0.992                                        10.5%
#         3                    0.933                                        14.9%
#        10                    0.688                                        30.8%
#        25                    0.529                                        37.3%
#
# At a loading of 1.0 the two indices correlate at 0.992: adjusting for the
# visible features very nearly controls the hidden one too, and the arm
# demonstrates nothing. Returns flatten past 10 because f0 is itself 56%
# predictable from f1-f11, capping how unobservable it can be made. 10 leaves
# roughly a third of targeting invisible to every estimator, which is the point
# of the arm.
HIDDEN_LOADING <- 10.0

# The targeting index is deliberately NONLINEAR in the features, so that a
# linear propensity model is genuinely misspecified rather than merely noisy:
#   g(X) = z + LAMBDA_SQ * (z^2 - 1) + KAPPA_INT * (x_a * x_b)
# where z is the standardized baseline-conversion log-odds index.
LAMBDA_SQ  <- 0.60
KAPPA_INT  <- 0.45
INTERACT_A <- "f3"
INTERACT_B <- "f8"

# ---- Observational subsample size ------------------------------------------
# The randomized benchmark is computed on all 13.98M rows. The confounded
# subsample is capped so that propensity matching and cross-fitted AIPW finish
# in minutes on a laptop. Override:  N_OBS=500000 Rscript run_all.R
N_OBS <- as.integer(Sys.getenv("N_OBS", "2000000"))
TREATED_SHARE_TARGET <- 0.60   # ad-exposed populations are the majority in practice

# ---- Cross-fitting ----------------------------------------------------------
N_FOLDS <- as.integer(Sys.getenv("N_FOLDS", "5"))

# ---- Estimator settings -----------------------------------------------------
PS_TRUNCATE  <- c(0.01, 0.99)  # propensity truncation for IPTW/AIPW
MATCH_CALIPER_SD <- 0.20       # caliper in SDs of the linear propensity score
N_BOOT_MATCH <- as.integer(Sys.getenv("N_BOOT_MATCH", "200"))

# ---- SLA guardrails ---------------------------------------------------------
# Declared here, before any estimate exists, and enforced by R/08_guardrails.R.
# An estimator that fails any of these is not fit to inform a budget decision,
# regardless of how close its point estimate looks.
GUARDRAILS <- list(
  # |estimate - benchmark| / |benchmark|
  max_relative_bias = 0.10,
  # the 95% CI must contain the randomized benchmark point estimate
  require_ci_covers_benchmark = TRUE,
  # observational CI must not be more than this multiple of the benchmark's
  max_ci_width_ratio = 3.0,
  # positivity: share of units outside PS_TRUNCATE must stay below this
  max_extreme_ps_share = 0.05,
  # covariate balance after adjustment, worst standardized mean difference
  max_abs_smd = 0.10,
  # effective sample size retained after weighting
  min_ess_ratio = 0.10
)

# ---- Helpers ----------------------------------------------------------------
db_connect <- function(read_only = FALSE) {
  # shared_home = FALSE only silences duckdb's extension-cache notice; this
  # pipeline loads no duckdb extensions.
  con <- DBI::dbConnect(duckdb::duckdb(dbdir = PATH_DUCKDB, read_only = read_only,
                                       shared_home = FALSE))
  DBI::dbExecute(con, "SET memory_limit='6GB'")
  DBI::dbExecute(con, "SET threads=6")
  con
}

read_sql <- function(path) paste(readLines(path, warn = FALSE), collapse = "\n")

write_table <- function(x, name) {
  dir.create(DIR_TABLES, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(x, file.path(DIR_TABLES, paste0(name, ".csv")), row.names = FALSE)
  invisible(x)
}
