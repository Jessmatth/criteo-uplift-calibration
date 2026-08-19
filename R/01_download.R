# ---- Fetch the Criteo uplift file -------------------------------------------
# The dataset is CC-BY-NC-SA 4.0. It is downloaded here rather than mirrored in
# the repository, both to respect the licence and to keep the checkout small.
source("R/00_config.R")

dir.create(DIR_RAW, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_DERIVED, showWarnings = FALSE, recursive = TRUE)

file_md5 <- function(p) as.character(tools::md5sum(p))

if (file.exists(PATH_GZ) && file.size(PATH_GZ) == DATA_BYTES) {
  message("Archive already present, skipping download.")
} else {
  message("Downloading ", DATA_URL, " (~297 MB) ...")
  # Criteo answers HEAD with 404 but honours GET; -L follows the 302 to Azure.
  ok <- utils::download.file(DATA_URL, PATH_GZ, mode = "wb", quiet = FALSE)
  if (ok != 0) stop("download failed")
}

stopifnot(file.size(PATH_GZ) == DATA_BYTES)
md5 <- file_md5(PATH_GZ)
if (!identical(md5, DATA_MD5)) {
  stop("MD5 mismatch: expected ", DATA_MD5, " got ", md5,
       "\nThe upstream file has changed; re-verify the row count before trusting results.")
}
message("Verified: ", DATA_BYTES, " bytes, MD5 ", md5)

utils::write.csv(
  data.frame(
    file       = basename(PATH_GZ),
    url        = DATA_URL,
    bytes      = DATA_BYTES,
    md5        = md5,
    rows       = DATA_ROWS,
    licence    = "CC-BY-NC-SA 4.0",
    citation   = "Diemert, Betlei, Renaudin, Amini (2018) A Large Scale Benchmark for Uplift Modeling. AdKDD@KDD.",
    retrieved  = format(Sys.Date())
  ),
  file.path(DIR_RAW, "MANIFEST.csv"), row.names = FALSE
)
message("Wrote ", file.path(DIR_RAW, "MANIFEST.csv"))
