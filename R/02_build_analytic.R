# ---- Stage into DuckDB and validate -----------------------------------------
source("R/00_config.R")
library(DBI)

con <- db_connect()

have <- "criteo" %in% dbListTables(con)
n_have <- if (have) dbGetQuery(con, "SELECT count(*) n FROM criteo")$n else 0L

if (!identical(as.numeric(n_have), as.numeric(DATA_ROWS))) {
  message("Staging ", basename(PATH_GZ), " into DuckDB ...")
  dbExecute(con, sprintf("SET VARIABLE gz_path = '%s'", PATH_GZ))
  invisible(dbExecute(con, read_sql("sql/01_stage.sql")))
} else {
  message("Table `criteo` already staged (", format(n_have, big.mark = ","), " rows).")
}

# ---- Validation -------------------------------------------------------------
# These are assertions about the experiment, not tests of the analysis. If any
# fail the upstream file is not what this protocol was written against.
n <- dbGetQuery(con, "SELECT count(*) n FROM criteo")$n
stopifnot(n == DATA_ROWS)

bad_exposure <- dbGetQuery(con,
  "SELECT count(*) n FROM criteo WHERE exposure = 1 AND treatment = 0")$n
stopifnot(bad_exposure == 0)   # exposure is reachable only under treatment

nulls <- dbGetQuery(con, sprintf(
  "SELECT %s FROM criteo",
  paste(sprintf("sum(CASE WHEN %s IS NULL THEN 1 ELSE 0 END) AS %s",
                c(FEATURES_ALL, "treatment", "conversion", "visit", "exposure"),
                c(FEATURES_ALL, "treatment", "conversion", "visit", "exposure")),
        collapse = ", ")))
stopifnot(all(unlist(nulls) == 0))

# ---- Profile ----------------------------------------------------------------
arm <- dbGetQuery(con, "
  SELECT treatment,
         count(*)          AS n,
         avg(conversion)   AS conversion_rate,
         avg(visit)        AS visit_rate,
         avg(exposure)     AS exposure_rate
  FROM criteo GROUP BY treatment ORDER BY treatment")
arm$arm <- ifelse(arm$treatment == 1, "treated (ad-eligible)", "control (held out)")
write_table(arm[c("arm", "n", "conversion_rate", "visit_rate", "exposure_rate")],
            "arm_summary")

cells <- dbGetQuery(con, "
  SELECT treatment, exposure, count(*) AS n,
         avg(conversion) AS conversion_rate, avg(visit) AS visit_rate
  FROM criteo GROUP BY 1,2 ORDER BY 1,2")
write_table(cells, "treatment_exposure_cells")

prof <- dbGetQuery(con, paste0("SELECT ", paste(sprintf(
  "min(%1$s) AS min_%1$s, max(%1$s) AS max_%1$s, avg(%1$s) AS mean_%1$s,
   stddev_samp(%1$s) AS sd_%1$s, count(DISTINCT %1$s) AS ndistinct_%1$s",
  FEATURES_ALL), collapse = ", "), " FROM criteo"))
profile <- data.frame(
  feature   = FEATURES_ALL,
  min       = as.numeric(prof[paste0("min_", FEATURES_ALL)]),
  max       = as.numeric(prof[paste0("max_", FEATURES_ALL)]),
  mean      = as.numeric(prof[paste0("mean_", FEATURES_ALL)]),
  sd        = as.numeric(prof[paste0("sd_", FEATURES_ALL)]),
  n_distinct = as.numeric(prof[paste0("ndistinct_", FEATURES_ALL)]),
  role = ifelse(FEATURES_ALL == FEATURE_HIDDEN, "hidden confounder", "visible to estimators")
)
write_table(profile, "feature_profile")

# Criteo coarsened the features, so distinct users collide on the full vector.
# This matters for matching: exact ties are common and tie-breaking must be
# explicit rather than left to sort order.
dupes <- dbGetQuery(con, "
  SELECT count(*) - count(DISTINCT (f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11)) AS duplicate_feature_rows
  FROM criteo")
message(sprintf("Rows sharing a feature vector with another row: %s (%.1f%%)",
                format(dupes$duplicate_feature_rows, big.mark = ","),
                100 * dupes$duplicate_feature_rows / n))

message("Staged and validated ", format(n, big.mark = ","), " rows.")

# Explicit rather than on.exit: these scripts are source()d in sequence by
# run_all.R, where an on.exit handler fires on source()'s frame and leaves the
# next script holding an invalid connection.
dbDisconnect(con, shutdown = TRUE)
