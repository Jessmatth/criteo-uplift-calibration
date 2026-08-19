# ---- SLA evaluation ---------------------------------------------------------
# Pure functions, unit-tested in tests/testthat/test-guardrails.R.
#
# The premise: an estimate is not usable because it is close to the truth. It is
# usable because it is close AND its uncertainty is honest AND the data could
# support the question in the first place. An estimator that lands on the
# benchmark with a 95% interval three times too wide, or by leaning on 2% of the
# sample, has not measured anything - it has guessed well once.

check_relative_bias <- function(estimate, benchmark, max_rel) {
  rel <- abs(estimate - benchmark) / abs(benchmark)
  list(value = rel, pass = is.finite(rel) && rel <= max_rel)
}

check_ci_covers <- function(lower, upper, benchmark, required) {
  covers <- lower <= benchmark && benchmark <= upper
  list(value = covers, pass = !required || covers)
}

check_ci_width <- function(lower, upper, bench_lower, bench_upper, max_ratio) {
  ratio <- (upper - lower) / (bench_upper - bench_lower)
  list(value = ratio, pass = is.finite(ratio) && ratio <= max_ratio)
}

check_positivity <- function(extreme_share, max_share) {
  list(value = extreme_share,
       pass = is.finite(extreme_share) && extreme_share <= max_share)
}

check_balance <- function(max_smd, limit) {
  list(value = max_smd, pass = is.finite(max_smd) && abs(max_smd) <= limit)
}

check_ess <- function(ess_ratio, min_ratio) {
  list(value = ess_ratio, pass = is.finite(ess_ratio) && ess_ratio >= min_ratio)
}

# Assemble the scorecard for one estimate. Checks that do not apply to an
# estimator (a naive difference has no propensity model, so positivity, balance
# and ESS are undefined) are recorded as NA and excluded from the verdict rather
# than silently passed.
evaluate_row <- function(est, bench, diag, g) {
  rb <- check_relative_bias(est$estimate, bench$benchmark, g$max_relative_bias)
  cv <- check_ci_covers(est$lower, est$upper, bench$benchmark,
                        g$require_ci_covers_benchmark)
  cw <- check_ci_width(est$lower, est$upper, bench$benchmark_lower,
                       bench$benchmark_upper, g$max_ci_width_ratio)
  has_ps <- !is.null(diag) && nrow(diag) == 1
  po <- if (has_ps) check_positivity(diag$extreme_ps_share, g$max_extreme_ps_share)
        else list(value = NA_real_, pass = NA)
  ba <- if (has_ps) check_balance(diag$max_abs_smd_weighted, g$max_abs_smd)
        else list(value = NA_real_, pass = NA)
  es <- if (has_ps) check_ess(diag$ess_ratio, g$min_ess_ratio)
        else list(value = NA_real_, pass = NA)

  checks <- list(relative_bias = rb, ci_covers_benchmark = cv,
                 ci_width_ratio = cw, positivity = po, balance = ba, ess = es)
  passes <- vapply(checks, function(x) as.logical(x$pass), logical(1))
  data.frame(
    relative_bias = rb$value, relative_bias_pass = rb$pass,
    ci_covers = cv$value, ci_covers_pass = cv$pass,
    ci_width_ratio = cw$value, ci_width_pass = cw$pass,
    extreme_ps_share = po$value, positivity_pass = po$pass,
    max_abs_smd_weighted = ba$value, balance_pass = ba$pass,
    ess_ratio = es$value, ess_pass = es$pass,
    checks_applicable = sum(!is.na(passes)),
    checks_failed = sum(!passes, na.rm = TRUE),
    verdict = if (all(passes, na.rm = TRUE)) "PASS" else "FAIL")
}
