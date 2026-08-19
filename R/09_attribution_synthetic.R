# =============================================================================
#  APPENDIX - SYNTHETIC MULTI-TOUCH JOURNEYS.  NOT REAL DATA.
# =============================================================================
# Everything in this file rests on user journeys that DO NOT EXIST in the Criteo
# log. Criteo publishes one row per user with a single ad channel and no
# touchpoint sequence, so first-touch, position-based and time-decay attribution
# are not computable from it. They are computable here only because this script
# invents a journey for each converter from a generative process written below.
#
# What is real: which users converted, and which users were exposed to the ad.
# Those come from the log. What is invented: how many other touchpoints each
# user had, on which channels, and in what order.
#
# Consequently the credit shares this file reports for organic, email and social
# are artifacts of the assumed process and mean nothing. The one number that
# survives is the comparison between what a multi-touch rule credits to DISPLAY
# and what the experiment says display actually caused - because the display
# touchpoint is placed from the real exposure flag, and the incremental total
# comes from the real randomized contrast.
#
# No result in README.md's headline tables comes from this file.
source("R/00_config.R"); source("R/estimator.R")
library(DBI)

set.seed(SEED)
con <- db_connect(read_only = TRUE)

# Only converters need journeys: attribution rules distribute credit FOR a
# conversion. Taking every treated converter rather than a sample keeps the
# totals exact.
conv <- dbGetQuery(con, sprintf(
  "SELECT user_id, exposure, %s FROM criteo WHERE treatment = 1 AND conversion = 1",
  paste(FEATURES_VISIBLE, collapse = ", ")))
n <- nrow(conv)

# ---- The invented generative process ----------------------------------------
# Intent index drives touch count: high-intent users are busier shoppers. The
# display ad is placed late in the journey when present, because retargeting
# fires after a user has already shown interest - which is precisely the
# mechanism that makes last-touch flatter display.
coefs <- readRDS("output/confound.rds")$coef
intent <- as.numeric(scale(as.matrix(conv[, FEATURES_VISIBLE]) %*%
                           coefs[FEATURES_VISIBLE]))
lambda <- exp(0.4 + 0.35 * intent)
n_other <- 1L + stats::rpois(n, lambda)

CHANNELS <- c("organic_search", "email", "paid_social", "direct")
CH_PROB  <- c(0.42, 0.22, 0.24, 0.12)
WINDOW_DAYS   <- 30
DISPLAY_WINDOW <- 10   # retargeting fires late, so display lands in the last third
HALFLIFE_DAYS <- 7

n_disp  <- as.integer(conv$exposure)
n_touch <- n_other + n_disp

# One flat table rather than a per-user rbind: at 36k converters the loop-and-
# bind version is minutes of work for an appendix.
user    <- c(rep(seq_len(n), n_other), rep(seq_len(n), n_disp))
channel <- c(sample(CHANNELS, sum(n_other), replace = TRUE, prob = CH_PROB),
             rep("display_ad", sum(n_disp)))
days    <- c(stats::runif(sum(n_other), 0, WINDOW_DAYS),
             stats::runif(sum(n_disp),  0, DISPLAY_WINDOW))
o <- order(user, -days)                       # within a user, oldest touch first
journeys <- data.frame(user = user[o], channel = channel[o], days_before = days[o])

k   <- n_touch[journeys$user]                                   # touches for this user
pos <- sequence(n_touch)                                        # 1..k within user

# ---- Credit rules -----------------------------------------------------------
first_touch <- as.numeric(pos == 1)
last_touch  <- as.numeric(pos == k)
linear      <- 1 / k
# 40/20/40: ends take 40% each, the middle splits 20%. Degenerate journeys of
# one or two touches fall back to 100% and 50/50.
position_based <- ifelse(k == 1, 1,
                  ifelse(k == 2, 0.5,
                  ifelse(pos == 1 | pos == k, 0.4, 0.2 / (k - 2))))
raw_decay  <- 2^(-journeys$days_before / HALFLIFE_DAYS)
time_decay <- raw_decay / ave(raw_decay, journeys$user, FUN = sum)

credits <- data.frame(channel = journeys$channel, first_touch, last_touch,
                      linear, position_based, time_decay)
RULES <- c("first_touch", "last_touch", "linear", "position_based", "time_decay")
stopifnot(all(abs(colSums(credits[RULES]) - n) < 1e-6))   # each rule splits exactly n conversions
by_channel <- aggregate(credits[RULES], list(channel = credits$channel), sum)
write_table(by_channel, "synthetic_channel_credit")

# ---- The one comparison that means anything ---------------------------------
cc <- utils::read.csv(file.path(DIR_TABLES, "attribution_credit.csv"))
cc <- cc[cc$outcome == "conversion", ]
disp <- by_channel[by_channel$channel == "display_ad", ]
cmp <- data.frame(
  rule = RULES,
  display_credited_conversions = as.numeric(disp[RULES]),
  incremental_conversions = cc$incremental_audited,
  overstatement = as.numeric(disp[RULES]) / cc$incremental_audited)
write_table(cmp, "synthetic_rule_comparison")

message("\n=== APPENDIX: synthetic journeys - channel mix is invented ===")
message(sprintf("  %d treated converters, %.2f touchpoints each on average.\n",
                n, nrow(journeys) / n))
message("  Conversions credited to DISPLAY by each rule, against the experiment:")
for (i in seq_len(nrow(cmp))) {
  message(sprintf("    %-16s %8s credited   %.2fx the %s that were incremental",
                  cmp$rule[i], format(round(cmp$display_credited_conversions[i]), big.mark = ","),
                  cmp$overstatement[i], format(round(cmp$incremental_conversions[i]), big.mark = ",")))
}
message("\n  Every rule is a rule for splitting observed conversions. None of them")
message("  can tell you how many would have happened anyway - only the holdout can.")

# Explicit rather than on.exit: these scripts are source()d in sequence by
# run_all.R, where an on.exit handler fires on source()'s frame and leaves the
# next script holding an invalid connection.
dbDisconnect(con, shutdown = TRUE)
